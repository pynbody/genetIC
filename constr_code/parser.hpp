//
// A simple-minded parser that takes commands like
//
//   command a b c
//
// and calls the appropriate function with arguments a,b,c
//
// Andrew Pontzen 2013


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

template <typename Rtype>
class Dispatch
{
private:

    class DispatchError: public std::runtime_error {
    public:
        DispatchError(const char* x) : runtime_error(x) { }

    };

    struct Base {
        virtual ~Base() {}
    };


    template<class R, class... Args>
    struct Func : Base
    {
        std::function<R(Args...)> f;
    };

    typedef std::function<Rtype(const Base*, std::istream&) > callerfn;

    std::map<std::string,
             std::pair<std::shared_ptr<Base>, callerfn> > _map;

    static void consume_comments(std::istream & input_stream) {
        std::string s;
        input_stream >> s;
        if(s[0]=='#' || s[0]=='%') {
            while(!input_stream.eof()) {
                input_stream >> s;
            }
        } else {
            if (s.size()!=0)
                throw DispatchError("Too many arguments");
        }
    }


    template<typename... Args>
    static Rtype call_function(const std::function<Rtype()> &f,
                                       std::istream & input_stream) {
       consume_comments(input_stream);
       return f();
    }

    template<typename T1>
    static Rtype call_function(const std::function<Rtype(T1)> &f,
                               std::istream & input_stream) {
        T1 t1;

        if(input_stream.eof())
            throw DispatchError("Insufficient arugments (expected 1, got 0)");
        input_stream >> t1;
        consume_comments(input_stream);
        return f(t1);
    }

    template<typename T1, typename T2>
    static Rtype call_function(const std::function<Rtype(T1,T2)> &f,
                               std::istream & input_stream) {

        if(input_stream.eof())
           throw DispatchError("Insufficient arugments (expected 2, got 0)");
        T1 t1;
        if(input_stream.eof())
           throw DispatchError("Insufficient arugments (expected 2, got 1)");
        T2 t2;
        input_stream >> t1 >> t2;
        consume_comments(input_stream);
        return f(t1,t2);
    }

    template<typename... Args>
    static Rtype unpack_and_call_function(const Base *fobj,
                               std::istream & input_stream) {

        // Turn the function back into the correct type (with argument information)
        // and call it using arguments on the input_stream

       auto pfunc = dynamic_cast<const Func<Rtype, Args...>*>(fobj);

       if(pfunc)
           return call_function<Args...>(pfunc->f, input_stream);
       else
           throw DispatchError("Wrong type");
    }

public:
    Dispatch() {}

    template<typename... Args>
    void add_route(const std::string &name,
              const std::function<Rtype(Args...)> &function) {
        // Define a route that calls the given function

        auto pfunc = std::make_shared<Func<Rtype, Args...>>();
        auto pcaller = std::function<Rtype(const Base*, std::istream&) >(&unpack_and_call_function<Args...>);
        pfunc->f = function;

        _map.insert(std::make_pair(name, std::make_pair(pfunc, pcaller)));
    }


    Rtype run(std::istream & input_stream) {
        // Take a single line out of the input stream and try to execute it

        std::string s;
        input_stream >> s;

        if (s[0]=='#' || s[0]=='%') // comment
            return;

        decltype(_map.at(s).first.get()) function;
        decltype(_map.at(s).second) caller;

        // retrieve the function and something to call it with
        try {
            function = _map.at(s).first.get();
            caller = _map.at(s).second;
        } catch (std::out_of_range &e) {
            throw DispatchError("Unknown command");
        }

        // do it!
        return caller(function, input_stream);
    }

    void run_loop(std::istream & input_stream) {
        std::string str;
        std::stringstream ss;

        int line =0;
        while (std::getline(input_stream, str)) {
            if(str.size()==0) continue;
            ss.str(str);
            line+=1;
            try {
                run(ss);
            } catch(DispatchError &e) {
                std::cerr << "Error \"" << e.what() << "\" on line " << line << " (\"" << str << "\")" << std::endl;
            }
            ss.clear();
        }

    }

};

template <typename Ctype, typename Rtype>
class ClassDispatch: public Dispatch<Rtype>
{
private:
    Ctype c;
public:
    ClassDispatch(Ctype &c) : c(c){}

    template<typename... Args>
    void add_class_route(const std::string &name, Rtype (Ctype::*f)(Args...)) {
        // make a lambda that performs the call
        auto call = [this,f](Args... input_args...) {(c.*f)(input_args...);};
        // add it as the route
        this->add_route(name, std::function<Rtype(Args...)>(call));
    }
};



/* TEST CODE --

class Stateful {
    double v=0;
public:
    void printout() {
        std::cerr << v << std::endl;
    }

    void add(double v1, double v2) {
        v=v1+v2;
    }

    void accum(double v1) {
        v+=v1;
    }
};

int main() {

    Stateful x;
    ClassDispatch<Stateful, void> dp(x);

    dp.add_class_route("printout", &Stateful::printout);
    dp.add_class_route("add", &Stateful::add);
    dp.add_class_route("accum", &Stateful::accum);

    dp.run_loop(std::cin);



}

*/
