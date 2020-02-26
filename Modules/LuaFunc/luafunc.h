// LuaFunc: simplifies use of runtime-provided lua-scripted functions.
// luafunc.h - API
// (C) 2013-2018 Nicholas G. Davies

// Example usage:
//
// #include <iostream>
// #include "luafunc.h"
// using namespace LuaFunc;
//
// int main()
// {
//     LuaFunction<double> L_sq("x", "x*x");        // function returns double: args [1] and expression [2]
//     LuaFunction<double> L_pow("x,y -> x^y");     // function returns double: compact notation, arrow separating args & expression
//     LuaFunction<double> L_fact("x -->> a = 1  for b = 2,x do a = a * b end  return a"); // big arrow notation, which allows fully specifying the function body
//
//     std::cout << L_sq(2) << "\n";     // prints 4
//     std::cout << L_pow(2, 3) << "\n"; // prints 8
//     std::cout << L_fact(5) << "\n";   // prints 120
// }

#ifndef LUAFUNC_H
#define LUAFUNC_H

#include <memory>
#include <lua/lua.hpp>
#include <string>
#include <iostream>
#include <stdexcept>

namespace LuaFunc
{
    namespace Internal
    {
        // The Lua state used by LuaFunc.
        extern lua_State* L;
        void Init();
        void Exit();

        // Push: pushes one or several values onto the Lua stack.
        inline void Push(int v)         { lua_pushinteger(L, v); }
        inline void Push(bool v)        { lua_pushboolean(L, v); }
        inline void Push(double v)      { lua_pushnumber(L, v); }
        inline void Push(const char* v) { lua_pushstring(L, v); }
        inline void Push(std::string v) { lua_pushstring(L, v.c_str()); }

        template <typename T, typename... U>
        void Push(T next, U... remaining)
        {
            Push(next);
            Push(remaining...);
        }

        // To: converts a value on the Lua stack into a C++ type.
        template <typename T>
        inline T To(int index) { T().CompilationError("Unregistered Type"); }

        template<> inline int To<int>(int index)                  { return lua_tointeger(L, index); }
        template<> inline bool To<bool>(int index)                { return lua_toboolean(L, index); }
        template<> inline double To<double>(int index)            { return lua_tonumber(L, index); }
        template<> inline const char* To<const char*>(int index)  { return lua_tostring(L, index); }
        template<> inline std::string To<std::string>(int index)  { return lua_tostring(L, index); }

        // RegistryEntry: an entry on the Lua registry, which cleans itself up upon destruction.
        // First, push the desired value onto the Lua stack, then call Register().
        class RegistryEntry
        {
        public:
            // Create a registry entry.
            RegistryEntry() : ref(LUA_NOREF) { }

            // Delete the registry entry.
            ~RegistryEntry()
            {
                if (Internal::L)
                    luaL_unref(Internal::L, LUA_REGISTRYINDEX, ref);
            }

            // Register the value at the top of the stack. Pops the value from the stack.
            void Register()
            {
                ref = luaL_ref(Internal::L, LUA_REGISTRYINDEX);
            }

            // Put the registered value at the top of the stack.
            void Get()
            {
                lua_rawgeti(Internal::L, LUA_REGISTRYINDEX, ref);
            }

        private:
            int ref;
        };
    };

    // LuaFunction<ReturnValue>: wraps a Lua function.
    template <typename ReturnValue>
    class LuaFunction
    {
    public:
        // default constructor
        LuaFunction() : n_args_pushed(0)
        { }

        // args of the form "x,y"; expr of the form "x^y+0.5"
        LuaFunction(std::string args, std::string expr)
         : compact_form(args + "->" + expr), n_args_pushed(0)
        {
            Create("return function(" + args + ") return " + expr + " end");
        }

        // compact_notation of the form x,y->x^y+0.5 or of the form 3^7+9 (no arguments, no arrow)
        LuaFunction(std::string compact_notation)
         : compact_form(compact_notation), n_args_pushed(0)
        {
            std::string::size_type arrow = compact_notation.find("->");
            std::string::size_type bigarrow = compact_notation.find("-->>");

            if (bigarrow != std::string::npos)
                Create("return function(" + compact_notation.substr(0, bigarrow) + ") " + compact_notation.substr(arrow + 4) + " end");
            else if (arrow != std::string::npos)
                Create("return function(" + compact_notation.substr(0, arrow) + ") return " + compact_notation.substr(arrow + 2) + " end");
            else
                Create("return function() return " + compact_notation + " end");
        }

        // Call function with provided arguments
        template <typename... T>
        ReturnValue operator()(T... args) const
        {
            lua_checkstack(Internal::L, sizeof...(args) + 1);   // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            Internal::Push(args...);                            // get arguments on stack
            lua_call(Internal::L, sizeof...(args), 1);          // call function, enforce one return value
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Call function with given vector as arguments
        template <typename ForwardIterator>
        ReturnValue Call(ForwardIterator first, ForwardIterator last) const
        {
            lua_checkstack(Internal::L, last - first + 1);      // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            for (ForwardIterator i = first; i != last; ++i)     // get arguments on stack
                Internal::Push(*i);
            lua_call(Internal::L, last - first, 1);             // call function, enforce one return value
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Push argument for later use of Call()
        template <typename Arg>
        LuaFunction<ReturnValue>& operator << (Arg a)
        {
            Internal::Push(a);
            ++n_args_pushed;
            return *this;
        }

        // Call function with previously-pushed arguments
        ReturnValue Call()
        {
            lua_checkstack(Internal::L, n_args_pushed + 1);     // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            lua_insert(Internal::L, -n_args_pushed - 1);        // move Lua function to correct stack position (before arguments)
            lua_call(Internal::L, n_args_pushed, 1);            // call function, enforce one return value
            n_args_pushed = 0;                                  // reset argument counter
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Call function with no arguments
        ReturnValue operator()()
        {
            return Call();
        }

        // Call function
        int CallMultRet()
        {
            int top = lua_gettop(Internal::L) - n_args_pushed;  // get stack index underneath arguments already pushed
            lua_checkstack(Internal::L, n_args_pushed + 128);   // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            lua_insert(Internal::L, -n_args_pushed - 1);        // move Lua function to correct stack position (before arguments)
            lua_call(Internal::L, n_args_pushed, LUA_MULTRET);  // call function, allowing multiple return values
            n_args_pushed = 0;                                  // reset argument counter
            n_return_values = lua_gettop(Internal::L) - top;    // get number of return values pushed to stack
            for (int n = n_return_values; n > 1; --n)
                lua_insert(Internal::L, -n);                    // reverse order of return values
            return n_return_values;                             // return number of return values pushed to stack
        }

        // Pop return value post-CallMultRet
        template <typename Ret>
        LuaFunction<ReturnValue>& operator >> (Ret& r)
        {
            if (n_return_values <= 0)
                throw std::runtime_error("No return values to pop.");
            r = Internal::To<Ret>(-1);                          // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            --n_return_values;
            return *this;
        }

        // Pop all remaining return values
        void FlushRet()
        {
            lua_pop(Internal::L, n_return_values);
            n_return_values = 0;
        }

        // Provide compact form of function definition
        std::string CompactForm() const
        {
            return compact_form;
        }

    private:
        std::shared_ptr<Internal::RegistryEntry> registry_entry;
        std::string compact_form;
        int n_args_pushed;
        int n_return_values;

        void Create(std::string command)
        {
            Internal::Init();                                   // initialise Lua, if needed
            registry_entry.reset(new Internal::RegistryEntry);  // create registry entry
            luaL_dostring(Internal::L, command.c_str());        // push function to stack
            registry_entry->Register();                         // store function in registry
        }
    };
};

#endif  // LUAFUNC_H