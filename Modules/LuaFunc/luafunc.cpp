// LuaFunc: simplifies use of runtime-provided lua-scripted functions.
// luafunc.cpp - implementation.

#include "luafunc.h"
#include <cstdlib>

namespace LuaFunc
{
    namespace Internal
    {
        // The Lua state used by LuaFunc.
        lua_State* L = 0;

        void Init()
        {
            // Register exit handler
            static bool ae = false;
            if (!ae)
            {
                ae = true;
                std::atexit(Exit);
            }

            // Create Lua state
            if (L == 0)
            {
                L = luaL_newstate();
                luaL_openlibs(L);
                lua_checkstack(L, 128);
            }
        }

        // Destroy all Lua states
        void Exit()
        {
            lua_close(L);
            L = 0;
        }
    };
};
