// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#ifndef __GEDIM_IOEnum_H
#define __GEDIM_IOEnum_H

#include <iostream>

namespace Gedim
{
namespace io_enum
{
// ***************************************************************************
template <auto V> constexpr std::string_view enum_name()
{
#if defined(_MSC_VER)
    constexpr std::string_view f = __FUNCSIG__;
    constexpr auto start = f.find_last_of('<') + 1;
    constexpr auto end = f.find_last_of('>');
#else
    constexpr std::string_view f = __PRETTY_FUNCTION__;
    constexpr auto start = f.find("V = ") + 4;
    constexpr auto end = f.find_first_of(";]", start);
#endif

    auto name = f.substr(start, end - start);

    if (auto p = name.find_last_of(':'); p != std::string_view::npos)
        name.remove_prefix(p + 1);

    return name;
}
// ***************************************************************************
template <auto V> constexpr bool is_valid()
{
    constexpr auto n = enum_name<V>();

    if (n.empty())
        return false;

    constexpr char c = n.front();

    return c == '_' || ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z');
}
// ***************************************************************************
template <typename E> struct enum_entry
{
    E value;
    std::string_view name;
};
// ***************************************************************************
template <typename E, int Min, std::size_t... I> constexpr auto make_table(std::index_sequence<I...>)
{
    constexpr auto count = (0 + ... + (is_valid<static_cast<E>(Min + int(I))>() ? 1 : 0));

    std::array<enum_entry<E>, count> table{};

    std::size_t j = 0;

    (
        [&] {
            constexpr auto v = static_cast<E>(Min + int(I));

            if constexpr (is_valid<v>())
                table[j++] = {v, enum_name<v>()};
        }(),
        ...);

    return table;
}
// ***************************************************************************
template <typename E, int Min = 0, int Max = 64> constexpr std::string_view enum_to_string(E value)
{
    constexpr auto table = make_table<E, Min>(std::make_index_sequence<Max - Min>{});

    for (auto const &e : table)
        if (e.value == value)
            return e.name;

    return "<unknown>";
}
// ***************************************************************************
} // namespace io_enum
} // namespace Gedim

#endif
