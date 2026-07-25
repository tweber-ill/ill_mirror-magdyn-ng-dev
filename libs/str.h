/**
 * string library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2013 - 2021
 * @note Forked on 7-Nov-2018 from my "tlibs" project (https://github.com/t-weber/tlibs).
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * Magpie
 * Copyright (C) 2022-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * TAS-Paths
 * Copyright (C) 2021       Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_STRINGS__
#define __TLIBS2_STRINGS__

#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <locale>
#include <limits>
#include <vector>
#include <stack>
#include <map>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <type_traits>
#include <cctype>
#include <cwctype>
#include <filesystem>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "expr.h"


namespace tl2 {
namespace fs = std::filesystem;


template<class t_str/*=std::string*/, class t_val/*=double*/>
std::pair<bool, t_val> eval_expr(const t_str& str) noexcept;


// -----------------------------------------------------------------------------


static inline std::wstring str_to_wstr(const std::string& str)
{
	return std::wstring(str.begin(), str.end());
}


static inline std::string wstr_to_str(const std::wstring& str)
{
	return std::string(str.begin(), str.end());
}


// overloaded in case the string is already of correct type
static inline const std::wstring& str_to_wstr(const std::wstring& str) { return str; }
static inline const std::string& wstr_to_str(const std::string& str) { return str; }


// -----------------------------------------------------------------------------



/**
 * removes all occurrences of a char in a string
 */
template<class t_str = std::string>
t_str remove_char(const t_str& str, typename t_str::value_type ch)
{
	t_str strRet;

	for(typename t_str::value_type c : str)
		if(c != ch)
			strRet.push_back(c);

	return strRet;
}


// -----------------------------------------------------------------------------


template<class t_str = std::string>
std::pair<t_str, t_str>
split_first(const t_str& str, const t_str& strSep, bool bTrim = false, bool bSeq = false)
{
	t_str str1, str2;

	std::size_t iLenTok = bSeq ? strSep.length() : 1;
	std::size_t ipos = bSeq ? str.find(strSep) : str.find_first_of(strSep);

	if(ipos != t_str::npos)
	{
		str1 = str.substr(0, ipos);
		if(ipos+iLenTok < str.length())
			str2 = str.substr(ipos+iLenTok, t_str::npos);
	}

	if(bTrim)
	{
		boost::trim(str1);
		boost::trim(str2);
	}

	return std::make_pair(str1, str2);
}


// ----------------------------------------------------------------------------


template<typename T, class t_str = std::string, bool bTIsStr = false>
struct _str_to_var_impl;


template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 1>
{
	inline const T& operator()(const t_str& str) const
	{
		return str;
	}
};


template<typename T, class t_str>
struct _str_to_var_impl<T, t_str, 0>
{
	inline T operator()(const t_str& str) const
	{
		if(!boost::trim_copy<t_str>(str).length())
			return T();

		T t{};
		typedef typename t_str::value_type t_char;
		std::basic_istringstream<t_char>{str} >> t;

		return t;
	}
};


/**
 * tokenises string on any of the chars in strDelim
 */
template<class T, class t_str = std::string, class t_cont = std::vector<T>>
void get_tokens(const t_str& str, const t_str& strDelim, t_cont& vecRet)
{
	using t_char = typename t_str::value_type;
	using t_tokeniser = boost::tokenizer<boost::char_separator<t_char>,
		typename t_str::const_iterator, t_str>;
	using t_tokiter = typename t_tokeniser::iterator;

	boost::char_separator<t_char> delim(strDelim.c_str());
	t_tokeniser tok(str, delim);

	for(t_tokiter iter = tok.begin(); iter != tok.end(); ++iter)
	{
		vecRet.push_back(
			_str_to_var_impl<T, t_str,
			std::is_convertible<T, t_str>::value>()(*iter));
	}
}


/**
 * tokenises string on strDelim
 */
template<class T, class t_str = std::string,
	template<class...> class t_cont = std::vector>
void get_tokens_seq(const t_str& str, const t_str& strDelim,
	t_cont<T>& vecRet, bool use_case = true)
{
	namespace algo = boost::algorithm;
	using t_char = typename t_str::value_type;

	std::locale loc;
	t_cont<t_str> vecStr;
	algo::iter_split(vecStr, str, algo::first_finder(strDelim,
		[use_case, &loc](t_char c1, t_char c2) -> bool
		{
			if(!use_case)
			{
				c1 = std::tolower(c1, loc);
				c2 = std::tolower(c2, loc);
			}

			return c1 == c2;
		}));

	for(const t_str& strTok : vecStr)
	{
		vecRet.push_back(
			_str_to_var_impl<T, t_str,
			std::is_convertible<T, t_str>::value>()(strTok));
	}
}


template<class T, class t_str = std::string, class t_cont = std::vector<T>>
bool parse_tokens(const t_str& str, const t_str& strDelim, t_cont& vecRet)
{
	std::vector<t_str> vecStrs;
	get_tokens<t_str, t_str, std::vector<t_str>>(str, strDelim, vecStrs);

	bool bOk = true;
	for(const t_str& str : vecStrs)
	{
		std::pair<bool, T> pairResult = eval_expr<t_str, T>(str);
		vecRet.push_back(pairResult.second);
		if(!pairResult.first)
			bOk = false;
	}

	return bOk;
}


template<typename T, class t_str = std::string>
T str_to_var_parse(const t_str& str)
{
	std::pair<bool, T> pairResult = eval_expr<t_str, T>(str);
	if(!pairResult.first)
		return T(0);
	return pairResult.second;
}


template<typename T, class t_str = std::string>
T str_to_var(const t_str& str)
{
	return _str_to_var_impl<T, t_str,
		std::is_convertible<T, t_str>::value>()(str);
}


/**
 * converts a string to a scalar value
 */
template<class t_scalar = double, class t_str = std::string>
t_scalar stoval(const t_str& str, bool pass_exception = false)
{
	try
	{
		if constexpr(std::is_same_v<t_scalar, float>)
			return std::stof(str);
		else if constexpr(std::is_same_v<t_scalar, double>)
			return std::stod(str);
		else if constexpr(std::is_same_v<t_scalar, long double>)
			return std::stold(str);
		else if constexpr(std::is_same_v<t_scalar, int>)
			return std::stoi(str);
		//else if constexpr(std::is_same_v<t_scalar, unsigned int>)
		//	return std::stoui(str);
		else if constexpr(std::is_same_v<t_scalar, long>)
			return std::stol(str);
		else if constexpr(std::is_same_v<t_scalar, unsigned long>)
			return std::stoul(str);
		else if constexpr(std::is_same_v<t_scalar, long long>)
			return std::stoll(str);
		else if constexpr(std::is_same_v<t_scalar, unsigned long long>)
			return std::stoull(str);

		// use the general conversion function if no specialised one was found
		return str_to_var<t_scalar, t_str>(str);
	}
	catch(const std::exception& ex)
	{
#ifdef __TLIBS2_SHOW_ERR__
		std::cerr << "String Lib: "
			<< "Conversion error: " << ex.what()
			<< "." << std::endl;
#endif
		if(pass_exception)
			throw ex;

		return t_scalar{};
	}
}


/**
 * unite bracket expressions in a token list to a single token
 */
template<class t_str = std::string, class t_cont = std::vector<t_str>>
t_cont unite_incomplete_tokens(const t_cont& toks, const t_str& opening = "([{", const t_str& closing = ")]}")
{
	using t_ch = typename t_str::value_type;

	std::stack<std::size_t> open_brackets;
	t_cont newtoks;

	// iterate tokens
	t_str cur_tok;
	for(auto tok_iter = toks.begin(); tok_iter != toks.end(); ++tok_iter)
	{
		// iterate token characters
		for(t_ch ch : *tok_iter)
		{
			bool has_opening_bracket = false;

			// iterate opening brackets
			for(std::size_t bracket_idx = 0; bracket_idx < opening.length(); ++bracket_idx)
			{
				t_ch bracket = opening[bracket_idx];

				if(ch == bracket)
				{
					open_brackets.push(bracket_idx);
					has_opening_bracket = true;
					break;
				}
			}

			if(has_opening_bracket)
				continue;

			// iterate closing brackets
			if(!open_brackets.empty())
			{
				std::size_t cur_bracket_idx = open_brackets.top();
				t_ch cur_bracket = closing[cur_bracket_idx];

				if(ch == cur_bracket)
					open_brackets.pop();
			}
		}

		if(cur_tok.length() != 0)
			cur_tok += " ";
		cur_tok += *tok_iter;

		if(open_brackets.empty())
		{
			newtoks.push_back(cur_tok);
			cur_tok = "";
		}
	}

	return newtoks;
}


// ----------------------------------------------------------------------------


template<class T, class t_ch,
	bool is_number_type = std::is_fundamental<T>::value>
struct _var_to_str_print_impl {};


template<class T, class t_ch> struct _var_to_str_print_impl<T, t_ch, false>
{
	void operator()(std::basic_ostream<t_ch>& ostr, const T& t) { ostr << t; }
};


template<class T, class t_ch> struct _var_to_str_print_impl<T, t_ch, true>
{
	void operator()(std::basic_ostream<t_ch>& ostr, const T& t)
	{
		// prevents printing "-0"
		T t0 = t;
		if(t0 == T(-0))
			t0 = T(0);

		ostr << t0;
	}
};


template<typename T, class t_str = std::string>
struct _var_to_str_impl
{
	t_str operator()(const T& t,
		std::streamsize iPrec = std::numeric_limits<T>::max_digits10,
		int iGroup = -1)
	{
		//if(std::is_convertible<T, t_str>::value)
		//	return *reinterpret_cast<const t_str*>(&t);

		typedef typename t_str::value_type t_char;

		std::basic_ostringstream<t_char> ostr;
		ostr.precision(iPrec);


		class Sep : public std::numpunct<t_char>
		{
		public:
			Sep() : std::numpunct<t_char>(1) {}
			~Sep() { /*std::cout << "~Sep();" << std::endl;*/ }
		protected:
			virtual t_char do_thousands_sep() const override { return ' ';}
			virtual std::string do_grouping() const override { return "\3"; }
		};
		Sep *pSep = nullptr;


		if(iGroup > 0)
		{
			pSep = new Sep();
			ostr.imbue(std::locale(ostr.getloc(), pSep));
		}

		_var_to_str_print_impl<T, t_char> pr;
		pr(ostr, t);
		t_str str = ostr.str();

		if(pSep)
		{
			ostr.imbue(std::locale());
			delete pSep;
		}
		return str;
	}
};


template<class t_str>
struct _var_to_str_impl<t_str, t_str>
{
	const t_str& operator()(const t_str& tstr, std::streamsize /*iPrec=10*/, int /*iGroup=-1*/)
	{
		return tstr;
	}

	t_str operator()(const typename t_str::value_type* pc, std::streamsize /*iPrec=10*/, int /*iGroup=-1*/)
	{
		return t_str(pc);
	}
};


template<typename T, class t_str=std::string>
t_str var_to_str(const T& t,
	std::streamsize iPrec = std::numeric_limits<T>::max_digits10-1,
	int iGroup = -1)
{
	_var_to_str_impl<T, t_str> _impl;
	return _impl(t, iPrec, iGroup);
}


/**
 * converts a container (e.g. a vector) to a string
 */
template<class t_cont, class t_str=std::string>
t_str cont_to_str(const t_cont& cont, const char* pcDelim=",",
	std::streamsize iPrec = std::numeric_limits<typename t_cont::value_type>::max_digits10-1)
{
	using t_elem = typename t_cont::value_type;

	t_str str;

	for(typename t_cont::const_iterator iter=cont.begin(); iter!=cont.end(); ++iter)
	{
		const t_elem& elem = *iter;

		str += var_to_str<t_elem, t_str>(elem, iPrec);
		if(iter+1 != cont.end())
			str += pcDelim;
	}
	return str;
}


// ----------------------------------------------------------------------------


template<class t_str /*=std::string*/, class t_val /*=double*/>
std::pair<bool, t_val> eval_expr(const t_str& str) noexcept
{
	if(boost::trim_copy<t_str>(str).length() == 0)
		return std::make_pair(true, t_val(0));

	try
	{
		ExprParser<t_val> parser;
		parser.SetAutoregisterVariables(false);
		bool ok = parser.parse(wstr_to_str(str));
		t_val valRes = parser.eval();
		return std::make_pair(ok, valRes);
	}
	catch(const std::exception& ex)
	{
#ifdef __TLIBS2_SHOW_ERR__
		std::cerr << "String Lib: "
			<< "Parsing failed with error: "
			<< ex.what() << "." << std::endl;
#endif
		return std::make_pair(false, t_val(0));
	}
}


// ----------------------------------------------------------------------------


/**
 * get the rgb colour values from a string
 */
template<class t_val>
bool get_colour(const std::string& _col, t_val *rgb)
{
	std::string col = boost::trim_copy(_col);

	// use default colour
	if(col == "" || col == "auto")
		return false;

	std::istringstream istrcolour(col);

	// optional colour code prefix
	if(istrcolour.peek() == '#')
		istrcolour.get();

	std::size_t colour = 0;
	istrcolour >> std::hex >> colour;

	if constexpr(std::is_floating_point_v<t_val>)
	{
		// get the colour values as floats in the range [0, 1]
		rgb[0] = t_val((colour & 0xff0000) >> 16) / t_val(0xff);
		rgb[1] = t_val((colour & 0x00ff00) >> 8) / t_val(0xff);
		rgb[2] = t_val((colour & 0x0000ff) >> 0) / t_val(0xff);
	}
	else
	{
		// get the colour values as bytes in the range [0, 255]
		rgb[0] = t_val((colour & 0xff0000) >> 16);
		rgb[1] = t_val((colour & 0x00ff00) >> 8);
		rgb[2] = t_val((colour & 0x0000ff) >> 0);
	}

	return true;
}


// ----------------------------------------------------------------------------


/**
 * get the extension of a file (without the point)
 */
template<class t_str = std::string>
t_str get_fileext(const t_str& file_path, bool to_lower = false)
{
	fs::path path(file_path);
	t_str file_ext = path.extension().string();

	// remove the point
	if(file_ext.length() > 0)
		file_ext = file_ext.substr(1);

	// convert to lower case
	if(to_lower)
		boost::to_lower(file_ext);

	return file_ext;
}


/**
 * get the file name without the extension
 */
template<class t_str = std::string>
t_str get_file_noext(const t_str& file_path, bool to_lower = false)
{
	fs::path path(file_path);

	fs::path path_noext = path.parent_path();
	path_noext /= path.stem();
	t_str pathname_noext = path_noext.string();

	// convert to lower case
	if(to_lower)
		boost::to_lower(pathname_noext);

	return pathname_noext;
}


/**
 * get the file name without the path
 */
template<class t_str = std::string>
t_str get_file_nodir(const t_str& file_path, bool to_lower = false)
{
	fs::path path(file_path);
	t_str file_stem = path.filename().string();

	// convert to lower case
	if(to_lower)
		boost::to_lower(file_stem);

	return file_stem;
}


}
#endif
