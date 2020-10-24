//
// A simple-minded parser that takes commands like
//
//   command a b c
//
// and calls the appropriate function with arguments a,b,c
//
// Andrew Pontzen 2013

#ifndef _PARSER_HPP
#define _PARSER_HPP

#include <memory>
#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <functional>
#include <iostream>
#include <exception>
#include <utility>
#include <functional>
#include <cctype>
#include <algorithm>
#include <vector>

#include "logging.hpp"

namespace tools {
  /*! \class DispatchError
      \brief Wrapper for runtime_errors, which specifically identifies errors related to parsing the input file.
  */
  class DispatchError : public std::runtime_error {
  public:
    //! Constructor from a char string - just construsts a runtime error.
    DispatchError(const char *x) : runtime_error(x) {}

  };

  //! Skips over any comments, and throws errors if we end up with too many arguments in an string.
  void consume_comments(std::istream &input_stream) {
    std::string s;
    input_stream >> s;
    if (s[0] == '#' || s[0] == '%') {
      while (!input_stream.eof()) {
        input_stream >> s;
      }
    } else {
      if (s.size() != 0)
        throw DispatchError("Too many arguments");
    }
  }

  template<typename Rtype, typename... Args>
  Rtype call_function(const std::function<Rtype(Args...)> &f, std::istream &input_stream,
                      std::ostream *output_stream);

  /*! \brief Consumes any comments, and attempts to run the function that template deduction has selected.
      \param f - function with no arguments
      \param input_stream - paramter file to be run

      Only works with functions with no arguments. Functions with arguments ultimately call this after
      recursively stripping the arguments using the specialisation below. This means we only call this
      function if we expect to be at the end of the string "function a b c ...".

      Consequently, the consume_comments sub-routine can assume there will be no more input arguments. If
      it encounters any text (ignoring them if they are comments), then it will throw an error. Otherwise,
      it will call the function.
  */
  template<typename Rtype>
  Rtype call_function(const std::function<Rtype()> &f,
                      std::istream &input_stream,
                      std::ostream *output_stream) {
    consume_comments(input_stream);
    if (output_stream) (*output_stream) << std::endl;
    return f();
  }

  /*! \brief Extracts an argument of the user specified function from the stream, and recurses the extract the other arguments
      \param f - function with first argument type T1, and an arbitrary list of other arguments
      \param input_stream - paramter file to be run
      \param output_stream - pointer to the stream used for output

      The function works by recursive template deduction. Given an input function with some arguments f(a,b,c,...),
      it constructs a new function, bound_f(b,c,...) which automatically calls f(a,b,c,...). This defines a new
      specialisation of the template for a function with one fewer arguments, and the operation is then defined
      recursively, extracting arguments one by one until we are left with bound_f(), which then triggers
      the specialisation above for functions with no arguments, ending the chain.
  */
  template<typename Rtype, typename T1, typename... Args>
  Rtype call_function(const std::function<Rtype(T1, Args...)> &f, std::istream &input_stream,
                      std::ostream *output_stream) {
    T1 arg1{};
    if (input_stream.eof()) {
      throw DispatchError("Insufficient arguments");
    } else {
      input_stream >> arg1;
      if (input_stream.fail())
        throw DispatchError("Error while parsing arguments");
    }

    if (output_stream != nullptr)
      (*output_stream) << arg1 << " ";

    // Create function with one fewer arguments:
    std::function<Rtype(Args...)> bound_f = [&f, arg1](auto &&... args) { return f(arg1, args...); };

    // Call the new function. Compiler will have to create a new specialisation of call_function:
    call_function<Rtype, Args...>(bound_f, input_stream, output_stream);

  };


  /*! \class Dispatch
      \brief Class that handles converting text strings like "function a b c" into function pointers that can be called, f(a,b,c)
  */
  template<typename Rtype>
  class Dispatch {
  private:

    /*! \struct Base
        \brief Base class of Func (allows us to pass Funcs with different argument lists to the same function)
    */
    struct Base {
      virtual ~Base() {}
    };


    /*! \struct Func
        \brief Struct to hold a single function, with arbitrary argument list and return type R
    */
    template<class R, class... Args>
    struct Func : Base {
      std::function<R(Args...)> f;
    };

    typedef std::function<Rtype(const Base *, std::istream &,
                                std::ostream *)> callerfn; //!< Wraps around a Func object, carrying additional information about the input and output stream

    std::map<std::string,
      std::pair<std::shared_ptr<Base>, callerfn> > _map; //!< Dictionary, used to map strings in the parameter file to functions that can be called.


    //! Turn the function back into the correct type (with argument information), and call it using arguments on the input_stream
    template<typename... Args>
    static Rtype unpack_and_call_function(const Base *fobj,
                                          std::istream &input_stream,
                                          std::ostream *output_stream) {

      auto pfunc = dynamic_cast<const Func<Rtype, Args...> *>(fobj);

      if (pfunc)
        return call_function<Rtype, Args...>(pfunc->f, input_stream, output_stream);
      else
        throw DispatchError("Wrong type");
    }

    //! Take a single line out of the input stream and try to execute it
    Rtype run(std::istream &input_stream, std::ostream *output_stream) {


      std::string s;
      input_stream >> s;

      if (s[0] == '#' || s[0] == '%') // comment
        return;

      decltype(_map.at(s).first.get()) function;
      decltype(_map.at(s).second) caller;

      std::transform(s.begin(), s.end(), s.begin(), ::tolower);

      // retrieve the function and something to call it with
      try {
        function = _map.at(s).first.get();
        caller = _map.at(s).second;
      } catch (std::out_of_range &e) {
        if (ignoreUnknownCommand)
          return Rtype();
        else
          throw DispatchError("Unknown command");
      }

      if (output_stream != nullptr)
        (*output_stream) << s << " ";

      // do it!
      return caller(function, input_stream, output_stream);
    }


    //! Parse the parameter file line by line, trying to run each command in sequence.
    void run_loop(std::istream &input_stream, std::ostream *output_stream) {
      std::string str;
      std::stringstream ss;

      int line = 0;
      while (std::getline(input_stream, str)) {
        line += 1;
        if (str.size() == 0) continue;
        ss.str(str);
        try {
          run(ss, output_stream);
        } catch (std::runtime_error &e) {
          logging::entry() << "Error \"" << e.what() << "\" on line " << line << " (\"" << str << "\")" << std::endl;
          exit(1);
        }
        ss.clear();
      }

    }

  public:
    bool ignoreUnknownCommand; //!< If true, the code won't throw an error if it doesn't recognise a command, and will just skip that line.

    //! Create a dispatcher, with ignoreUnknownCommand specified.
    Dispatch(bool ignoreUnknownCommand = false) : ignoreUnknownCommand(ignoreUnknownCommand) {}

    /*! \brief Define a route that calls the given function
        \param name - string that will be associated to the specified function
        \param function - function that will be called for a given string

        This function essentially ties together a string and a function pointer, and adds them to
        the _map object stored in the Dispatch class. When the parser attempts to run functions,
        it will look through this map to see if any match.

        Note that there is a different add_route function for each pair of string and
        function pointer, created by the compiler. To add new commands, these
        functions need to be called in main.cpp, so that the compiler knows to
        create them and the map can be created at run-time.
    */
    template<typename... Args>
    void add_route(const std::string &name,
                   const std::function<Rtype(Args...)> &function) {

      auto pfunc = std::make_shared<Func<Rtype, Args...>>();
      auto pcaller = std::function<Rtype(const Base *, std::istream &, std::ostream *)>(
        &unpack_and_call_function<Args...>);
      pfunc->f = function;

      std::string lname(name);

      std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);
      _map.insert(std::make_pair(lname, std::make_pair(pfunc, pcaller)));
    }

    //! Converts an input string to an istringstream and attempts to run it
    Rtype run(std::string input) {
      std::istringstream ss(input);
      run(ss);
    }

    //! Runs an input stream with no output stream
    void run(std::istream &input_stream) {
      run(input_stream, nullptr);
    }

    //! Runs an input stream with an output stream specified
    void run(std::istream &input_stream, std::ostream &output_stream) {
      run(input_stream, &output_stream);
    }

    //! Run a parsing loop over the parameter file with the specified input stream, and no output stream
    void run_loop(std::istream &input_stream) {
      run_loop(input_stream, nullptr);
    }

    //! Run a parsing loop over the parameter file with the specified input and output stream
    void run_loop(std::istream &input_stream, std::ostream &output_stream) {
      run_loop(input_stream, &output_stream);
    }

  };

  template<typename Ctype, typename Rtype>
  class InstanceDispatch;

  /*! \class ClassDispatch
      \brief Allows us to map class member function into a pair with a string, so that these can be run by the parser.

      This is necessary, because member functions are technically of a
      different type to ordinary functions (they have an
      implicit 'this' parameter, for example). We need to bind them with
      the instance of the class that will actually call them and convert
      the result into std::function objects that the add_route member
      function of the Dispatch class can understand.

      Note that we cannot simply setup all these functions when
      add_class_route is first called in main.cpp, because member functions
      can only be called if we have an instance of the class to call
      them with (which in the case of using an input mapper, may not
      have been created yet). Thus we split the task in two:

      add_class_route - adds to the list of (string,function) pairs
                        that are needed for a given instance
                        (stored in the adderFunctions vector)
      specify_instance - for a specified instance, creates the
                         std::functions that can then be added to
                         the (string,function) map and called.

      What is then stored in the ClassDispatch object is actually a
      list of functions which take an InstanceDispatch object defined for
      the required instance as their argument. They then bind it to the
      appropriate member function, creating the (string,function) pair,
      and add it to the list of callable functions.
  */
  template<typename Ctype, typename Rtype>
  class ClassDispatch {
  private:
    std::vector<std::function<void(
      InstanceDispatch<Ctype, Rtype> &)>> adderFunctions; //!< List of functions that need to be added for a given instance of the class
  public:
    //! Default constructor
    ClassDispatch() {}

    //! Adds a (string,function) pair to the list of pairs that need to be setup for a given instance of the Ctype class.
    template<typename... Args>
    void add_class_route(const std::string &name, Rtype (Ctype::*f)(Args...)) {
      auto addcall = std::function<void(InstanceDispatch<Ctype, Rtype> &)>(
        [name, f](InstanceDispatch<Ctype, Rtype> &pDispatchObj) {
          pDispatchObj.add_class_route(name, f);
        });


      // make a lambda that adds the route to a specific object
      adderFunctions.emplace_back(addcall);
    }

    //! Adds a (string,function) pair to the list of pairs that need to be setup for a given instance of the Ctype class.
    template<typename... Args>
    void add_deprecated_class_route(const std::string &name, const std::string &preferredName, Rtype (Ctype::*f)(Args...)) {
      auto addcall = std::function<void(InstanceDispatch<Ctype, Rtype> &)>(
        [name, preferredName, f](InstanceDispatch<Ctype, Rtype> &pDispatchObj) {
          pDispatchObj.add_deprecated_class_route(name, preferredName, f);
        });

      adderFunctions.emplace_back(addcall);
    }

    //! For the given instance of Ctype, setup class routes for its member functions
    InstanceDispatch<Ctype, Rtype> specify_instance(Ctype &c) {
      auto id = InstanceDispatch<Ctype, Rtype>(c);
      for (auto fn : adderFunctions)
        fn(id); // Calls the add_class_route member member function of InstanceDispatch, creating a (string,function) pair and adding it to the map
      return id;

    }
  };

  /*! \class InstanceDispatch
      \brief Derived class of Dispatch, specialised to be able to add member functions of another class as members of a (string,function) pair.
  */
  template<typename Ctype, typename Rtype>
  class InstanceDispatch : public Dispatch<Rtype> {
  private:
    Ctype *pC; //!< Pointer to the underlying instance of the relevant class
    ClassDispatch<Ctype, Rtype> prototype;

  public:
    //! Constructor from an instance of the specified class
    InstanceDispatch(Ctype &c) : pC(&c) {}

    //! Wraps the member function inside a std::function that can be added to the (string,function) map.
    template<typename... Args>
    void add_class_route(const std::string &name, Rtype (Ctype::*f)(Args...)) {
      // make a lambda that performs the call
      auto call = [this, f](Args... input_args) { (pC->*f)(input_args...); };
      // add it as the route
      this->add_route(name, std::function<Rtype(Args...)>(call));
    }

    //! As add_class_route, but emits a deprecation warning if the call is made
    template<typename... Args>
    void add_deprecated_class_route(const std::string &name, const std::string &preferredName, Rtype (Ctype::*f)(Args...)) {
      auto call = [this, f, name, preferredName](Args... input_args) {
        logging::entry(logging::level::warning) << "WARNING: " << name << " is a deprecated command and has been replaced by " << preferredName
                  << std::endl;
        (pC->*f)(input_args...);
      };

      this->add_route(name, std::function<Rtype(Args...)>(call));
    }


  };

}
#endif
