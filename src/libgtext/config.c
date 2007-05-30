/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/config.h>
#include <assert.h>
#include <string.h>
#include <libgtcore/log.h>
#include <libgtcore/ensure.h>
#include <libgtext/config.h>
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"

struct Config
{
	lua_State *L;
  Str *fn;
};

static void gtlua_new_table(lua_State *L, const char *key)
{
  lua_pushstring(L, key);
	lua_newtable(L);
	lua_settable(L, -3);
}

#define ADDCOLOR(rgb) \
	lua_pushstring(cfg->L, #rgb); \
	lua_pushnumber(cfg->L, color.rgb); \
	lua_settable(cfg->L, -3);

#define GET_AND_SET_COLOR_VALUE(rgb) \
  lua_getfield(cfg->L, -1, #rgb); \
  if (lua_isnil(cfg->L, -1) || !lua_isnumber(cfg->L, -1)) { \
    env_log_log(env, "%s color value for type '%s' is undefined " \
		"or not numeric, using default",#rgb, key); \
  } \
  else \
  { \
    color.rgb = lua_tonumber(cfg->L,-1); \
  } \
  lua_pop(cfg->L, 1); \

/*!
Creates a Config object.
\param env Pointer to Environment object.
\return Pointer to the new object.
*/
Config* config_new(Env *env)
{
	Config *cfg;
  env_error_check(env);
	/* open and assign log facility */
	Log* log = log_new(env_ma(env));
	env_set_log(env, log);
  cfg = env_ma_malloc(env, sizeof (Config));
	cfg->fn = NULL;
  cfg->L = luaL_newstate();
  if (!cfg->L)
  {
    env_error_set(env, "out of memory (cannot create new lua state)");
  }
	else luaL_openlibs(cfg->L);
  return cfg;
}

/*!
Deletes a Config object.
\param cfg Pointer to Config object to delete.
\param env Pointer to Environment object.
*/
void config_delete(Config *cfg, Env *env)
{
  assert(cfg);
	if (cfg->L) lua_close(cfg->L);
	if (cfg->fn) str_delete(cfg->fn,env);
  env_ma_free(cfg, env);
}

/*!
Loads and executes a Lua configuration file.
This file must contain a table called 'config'.
\param cfg Config object to load into.
\param fn Filename of the script to execute.
\param env Pointer to Environment object.
*/
void config_load_file(Config *cfg, Str *fn, Env* env)
{
  int has_err = 0;
  env_error_check(env);
  assert(cfg && cfg->L && fn);
  cfg->fn = str_ref(fn);
  if (luaL_loadfile(cfg->L, str_get(fn)) ||
      lua_pcall(cfg->L, 0, 0, 0))
	{
    env_error_set(env, "cannot run configuration file: %s",
                  lua_tostring(cfg->L, -1));
    has_err = -1;
  }
  if (!has_err)
	{
    lua_getglobal(cfg->L, "config");
    if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1))
		{
      env_error_set(env, "'config' is not defined or not a table in \"%s\"",
                str_get(fn));
    }
		lua_pop(cfg->L, 1);
  }
}

/*!
Searches for a section inside the config table,
creating it if it does not exist and finally pushing
it on the Lua stack (at the top).
\param cfg Config object to search in.
\param section Name of the section to search for.
\param env Pointer to Environment object.
\return Total number of items pushed on the stack by this function.
*/
static int config_find_section_for_setting(Config* cfg,
                                           const char* section,
																					 Env* env)
{
  int depth = 0;
  assert(cfg && section);
	env_error_check(env);
	lua_getglobal(cfg->L, "config");
	if (lua_isnil(cfg->L, -1))
	{
	  lua_newtable(cfg->L);
		lua_setglobal(cfg->L, "config");
		lua_getglobal(cfg->L, "config");
	}
	depth++;
	lua_getfield(cfg->L, -1, section);
	if (lua_isnil(cfg->L, -1))
	{
	  lua_pop(cfg->L, 1);
    gtlua_new_table(cfg->L, section);
		lua_getfield(cfg->L, -1, section);
	}
	depth++;
	return depth;
}

/*!
Searches for a section inside the config table,
returning if it is not found.
\param cfg Config object to search in.
\param section Name of the section to search for.
\param env Pointer to Environment object.
\return -1 if not found, otherwise number of items
        pushed on the stack by this function.
*/
static int config_find_section_for_getting(Config* cfg,
                                           const char* section,
																					 Env* env)
{
  int depth = 0;
	assert(cfg && section);
	lua_getglobal(cfg->L, "config");
	if (lua_isnil(cfg->L, -1))
	{
	  env_log_log(env, "'config' is not defined");
		lua_pop(cfg->L, 1);
	  return -1;
	} else depth++;
	lua_getfield(cfg->L, -1, section);
	if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1))
	{
    env_log_log(env, "section '%s' is not defined", section);
		lua_pop(cfg->L, 1);
		return -1;
	} else depth++;
	return depth;
}

/*!
Sets a string value in the configuration to a certain value.
\param cfg Config object to search in.
\param section Section to set a key in.
\param key Key to set a value for.
\param str String that is to be set.
\param env Pointer to Environment object.
*/
void config_set_cstr(Config *cfg,
                     const char* section,
										 const char *key,
										 const char* str,
										 Env* env)
{
  int i = 0;
	assert(cfg && section && key);
  i = config_find_section_for_setting(cfg, section, env);
	lua_pushstring(cfg->L, key);
	lua_pushstring(cfg->L, str);
	lua_settable(cfg->L, -3);
	lua_pop(cfg->L, i);
}

/*!
Sets a color value in the configuration to a certain value.
\param cfg Config object to search in.
\param key Key (e.g. feature) to set a color for.
\param color Color to associate with key.
\param env Pointer to Environment object.
*/
void config_set_color(Config *cfg, const char *key, Color color, Env* env)
{
  int i = 0;
	assert(cfg && key);
  i = config_find_section_for_setting(cfg, "colors", env);
	lua_getfield(cfg->L, -1, key);
	i++;
	if (lua_isnil(cfg->L, -1))
	{
	  lua_pop(cfg->L, 1);
    gtlua_new_table(cfg->L, key);
		lua_getfield(cfg->L, -1, key);
	}
  ADDCOLOR(red);
	ADDCOLOR(green);
	ADDCOLOR(blue);
	lua_pop(cfg->L, i);
}

/*!
Retrieves a string value from the configuration.
\param cfg Config object to search in.
\param section Section to get a key from.
\param key Key to get a value from.
\param env Pointer to Environment object.
\return NULL if not found, else string pointer.
*/
const char* config_get_cstr(Config *cfg,
                            const char* section,
														const char *key,
														Env* env)
{
  assert(cfg && key && section);
  const char* str = NULL;
	int i = 0;
 	env_error_check(env);
	/* get section */
	i = config_find_section_for_getting(cfg, section, env);
	/* could not get section, return default */
  if (i < 0) {
	  return str;
	}
	/* lookup entry for given key */
  lua_getfield(cfg->L, -1, key);
	if (lua_isnil(cfg->L, -1) || !lua_isstring(cfg->L, -1))
	{
    env_log_log(env, "no value is defined for key '%s'",
		            key);
		lua_pop(cfg->L, 1);
		return str;
  } else i++;
	/* retrieve string */
  str = lua_tostring(cfg->L, -1);
	/* reset stack to original state for subsequent calls */
	lua_pop(cfg->L, i);
	return str;
}

/*!
Retrieves a color value from the configuration.
\param cfg Config object to search in.
\param key Key (e.g. feature) to get a color for.
\param env Pointer to Environment object.
\return color Color associated with key.
*/
Color config_get_color(Config *cfg, const char *key, Env* env)
{
  assert(cfg && key);
  Color color;
	int i = 0;
 	env_error_check(env);
	/* set default colors */
	color.red=0.8; color.green = 0.8; color.blue=0.8;
	/* get section */
	i = config_find_section_for_getting(cfg, "colors", env);
	/* could not get section, return default */
  if (i < 0) return color;
	/* lookup color entry for given feature */
  lua_getfield(cfg->L, -1, key);
	if (lua_isnil(cfg->L, -1) || !lua_istable(cfg->L, -1))
	{
    env_log_log(env, "no colors are defined for type '%s', will use defaults",
		            key);
		lua_pop(cfg->L, 1);
		return color;
  } else i++;
	/* update color struct */
  GET_AND_SET_COLOR_VALUE(red);
	GET_AND_SET_COLOR_VALUE(green);
	GET_AND_SET_COLOR_VALUE(blue);
	/* reset stack to original state for subsequent calls */
	lua_pop(cfg->L, i);
	return color;
}

/*!
Reloads the Lua configuration file.
\param cfg Config object to search in.
\param env Pointer to Environment object.
*/
void config_reload(Config *cfg, Env *env)
{
  assert(cfg && cfg->fn);
  config_load_file(cfg, cfg->fn, env);
}

/*!
Unit tests for the Config class.
\param env Pointer to Environment object.
\return Error status.
*/
int config_unit_test(Env* env)
{
  int has_err = 0;
	Config *cfg;
	const char* test1 = "mRNA";
	const char* test2 = "gene";
  const char* str;
	Str *luafile = str_new_cstr("config.lua",env);
	Color col1, col2, col, defcol, tmpcol;

  /* example colors */
  col1.red=.1;col1.green=.2;col1.blue=.3;
  col2.red=.4;col2.green=.5;col2.blue=.6;
  col.red=1;col.green=1;col.blue=1;
	defcol.red=.8;defcol.green=.8;defcol.blue=.8;

  /* instantiate new config object */
	cfg = config_new(env);

	/* at the beginning, all values are defaults, since nothing is defined */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "cds", env);
  ensure(has_err, color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, color_equals(tmpcol,defcol));
  str = config_get_cstr(cfg, "collapse", "exon", env);
  ensure(has_err, (str == NULL));

  /* execute the config file */
	config_load_file(cfg, luafile, env);

  /* now we expect the values to exist and have certain values */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, !color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, color_equals(tmpcol,col1));
  tmpcol = config_get_color(cfg, "cds", env);
  ensure(has_err, color_equals(tmpcol,col2));
	tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, color_equals(tmpcol,defcol));
  str = config_get_cstr(cfg, "collapse", "exon", env);
  ensure(has_err, (str != NULL));
  ensure(has_err, (strcmp(str,test1)==0));
	str = config_get_cstr(cfg, "bar", "baz", env);
  ensure(has_err, (str == NULL));

  /* change some values... */
  config_set_color(cfg, "exon", col, env);
  config_set_cstr(cfg, "collapse", "exon", test2, env);

	/* is it saved correctly? */
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, !color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "exon", env);
  ensure(has_err, color_equals(tmpcol,col));
  str = config_get_cstr(cfg, "collapse", "exon", env);
  ensure(has_err, (str != NULL));
  ensure(has_err, !(strcmp(str,test1)==0));
  ensure(has_err, (strcmp(str,test2)==0));

  /* create a new color definition */
  config_set_color(cfg, "foo", col, env);
  config_set_cstr(cfg, "bar", "baz", test1, env);

	/* is it saved correctly? */
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, !color_equals(tmpcol,defcol));
  tmpcol = config_get_color(cfg, "foo", env);
  ensure(has_err, color_equals(tmpcol,col));
	str = config_get_cstr(cfg, "bar", "baz", env);
  ensure(has_err, (str != NULL));
  ensure(has_err, (strcmp(str,test1)==0));
	str = config_get_cstr(cfg, "bar", "test", env);
  ensure(has_err, (str == NULL));

  /* mem cleanup */
  str_delete(luafile, env);
	config_delete(cfg, env);

	return has_err;
}
