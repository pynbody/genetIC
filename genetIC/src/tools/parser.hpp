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

namespace tools {
  class DispatchError : public std::runtime_error {
  public:
    DispatchError(const char *x) : runtime_error(x) {}

  };

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
  Rtype call_function(const std::function<Rtype( Args...)> &f, std::istream &input_stream,
                             std::ostream * output_stream);

  template<typename Rtype>
  Rtype call_function(const std::function<Rtype()> &f,
                             std::istream &input_stream,
                             std::ostream * /* *output_stream*/) {
    consume_comments(input_stream);
    return f();
  }

  template<typename Rtype, typename T1, typename... Args>
  Rtype call_function(const std::function<Rtype(T1, Args...)> &f, std::istream &input_stream,
                             std::ostream * output_stream) {
    T1 arg1;
    if (input_stream.eof())
      throw DispatchError("Insufficient arugments");
    input_stream >> arg1;

    if (output_stream != nullptr)
      (*output_stream) << arg1 << " " << std::endl;

    std::function<Rtype(Args...)> bound_f = [&f, arg1](auto&&... args) { return f(arg1, args...); };

    call_function<Rtype, Args...>(bound_f, input_stream, output_stream);

  };


  template<typename Rtype>
  class Dispatch {
  private:

    struct Base {
      virtual ~Base() {}
    };


    template<class R, class... Args>
    struct Func : Base {
      std::function<R(Args...)> f;
    };

    typedef std::function<Rtype(const Base *, std::istream &, std::ostream *)> callerfn;

    std::map<std::string,
      std::pair<std::shared_ptr<Base>, callerfn> > _map;


    template<typename... Args>
    static Rtype unpack_and_call_function(const Base *fobj,
                                          std::istream &input_stream,
                                          std::ostream *output_stream) {

      // Turn the function back into the correct type (with argument information)
      // and call it using arguments on the input_stream

      auto pfunc = dynamic_cast<const Func<Rtype, Args...> *>(fobj);

      if (pfunc)
        return call_function<Rtype, Args...>(pfunc->f, input_stream, output_stream);
      else
        throw DispatchError("Wrong type");
    }

    Rtype run(std::istream &input_stream, std::ostream *output_stream) {
      // Take a single line out of the input stream and try to execute it

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
        } catch (std::exception &e) {
          std::cerr << "Error \"" << e.what() << "\" on line " << line << " (\"" << str << "\")" << std::endl;
          exit(1);
        }
        ss.clear();
      }

    }

  public:
    bool ignoreUnknownCommand;

    Dispatch(bool ignoreUnknownCommand = false) : ignoreUnknownCommand(ignoreUnknownCommand) {}

    template<typename... Args>
    void add_route(const std::string &name,
                   const std::function<Rtype(Args...)> &function) {
      // Define a route that calls the given function

      auto pfunc = std::make_shared<Func<Rtype, Args...>>();
      auto pcaller = std::function<Rtype(const Base *, std::istream &, std::ostream *)>(
        &unpack_and_call_function<Args...>);
      pfunc->f = function;

      std::string lname(name);

      std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);
      _map.insert(std::make_pair(lname, std::make_pair(pfunc, pcaller)));
    }

    Rtype run(std::string input) {
      std::istringstream ss(input);
      run(ss);
    }

    void run(std::istream &input_stream) {
      run(input_stream, nullptr);
    }

    void run(std::istream &input_stream, std::ostream &output_stream) {
      run(input_stream, &output_stream);
    }

    void run_loop(std::istream &input_stream) {
      run_loop(input_stream, nullptr);
    }

    void run_loop(std::istream &input_stream, std::ostream &output_stream) {
      run_loop(input_stream, &output_stream);
    }

  };

  template<typename Ctype, typename Rtype>
  class InstanceDispatch;

  template<typename Ctype, typename Rtype>
  class ClassDispatch {
  private:
    std::vector<std::function<void(InstanceDispatch<Ctype, Rtype> &)>> adderFunctions;
  public:
    ClassDispatch() {}

    template<typename... Args>
    void add_class_route(const std::string &name, Rtype (Ctype::*f)(Args...)) {
      auto addcall = std::function<void(InstanceDispatch<Ctype, Rtype> &)>(
        [name, f](InstanceDispatch<Ctype, Rtype> &pDispatchObj) {
          pDispatchObj.add_class_route(name, f);
        });


      // make a lambda that adds the route to a specific object
      adderFunctions.emplace_back(addcall);
    }

    InstanceDispatch<Ctype, Rtype> specify_instance(Ctype &c) {
      auto id = InstanceDispatch<Ctype, Rtype>(c);
      for (auto fn : adderFunctions)
        fn(id);
      return id;

    }
  };

  template<typename Ctype, typename Rtype>
  class InstanceDispatch : public Dispatch<Rtype> {
  private:
    Ctype *pC;
    ClassDispatch<Ctype, Rtype> prototype;

  public:
    InstanceDispatch(Ctype &c) : pC(&c) {}

    template<typename... Args>
    void add_class_route(const std::string &name, Rtype (Ctype::*f)(Args...)) {
      // make a lambda that performs the call
      auto call = [this, f](Args... input_args) { (pC->*f)(input_args...); };
      // add it as the route
      this->add_route(name, std::function<Rtype(Args...)>(call));
    }


  };

}
#endif