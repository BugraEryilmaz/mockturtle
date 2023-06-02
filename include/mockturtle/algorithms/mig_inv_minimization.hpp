/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file mig_inv_minimization.hpp
  \brief MIG Fanout Restriction Algorithm

  \author Bugra Eryilmaz
*/

#pragma once

#include "../networks/mig.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/fanout_view.hpp"

#include <iostream>
#include <optional>
#include <unordered_map>

#define DLEVEL 0
#if DLEVEL == 0
#define DEBUG_PRINT_MIG( x )
#define DEBUG( x )
#define INFO( x )
#elif DLEVEL == 1
#define DEBUG_PRINT_MIG( x )
#define DEBUG( x )
#define INFO( x ) std::cout << x << std::endl
#elif DLEVEL == 2
#include <mockturtle/utils/debugging_utils.hpp>
#define DEBUG_PRINT_MIG( x ) print( x )
#define DEBUG( x ) std::cout << x << std::endl
#define INFO( x ) std::cout << x << std::endl
#endif

namespace mockturtle
{

/*! \brief Parameters for mig_inv_minimization.
 *
 * The data structure `mig_inv_minimization_params` holds configurable
 * parameters with default arguments for `mig_inv_minimization`.
 */
struct mig_inv_minimization_params
{
};

/*! \brief Statistics for mig_inv_minimization. */
struct mig_inv_minimization_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Total number of recursive calls. */
  uint64_t num_calls{ 0 };

  /*! \brief Increase in the number of nodes. */
  uint64_t num_inverters_removed{ 0 };
};

namespace detail
{
template<class Ntk>
class mig_inv_minimization_impl
{
public:
  mig_inv_minimization_impl( Ntk& ntk, mig_inv_minimization_params const& ps, mig_inv_minimization_stats& st )
      : ntk( ntk ), ps( ps ), st( st )
  {
  }
  void run()
  {
    stopwatch t( st.time_total );
    bool minimized = true;
    while ( minimized )
    {
      minimized = false;
      ntk.foreach_gate( [this, &minimized]( auto const& n ) {
        INFO( "node: " << ntk.node_to_index( n ) );
        if ( ntk.is_dead( n ) )
        {
          INFO( "    dead node" );
          return;
        }
        if ( int gain = one_level( n ); gain > 0 )
        {
          INFO( "    one level inverted" );
          st.num_inverters_removed += gain;
          inv_node( n );
          DEBUG_PRINT_MIG( ntk );
          minimized = true;
        }
        if ( int gain = two_level( n ); gain > 0 )
        {
          INFO( "    two level inverted" );
          st.num_inverters_removed += gain;
          auto new_node = inv_node( n );
          ntk.foreach_fanout( new_node, [&, &minimized]( auto const& fo ) {
            INFO( "  Parent node: " << ntk.node_to_index( fo ) );
            if ( int parent_gain = one_level( fo ); parent_gain > 0 )
            {
              INFO( "    one level inverted" );
              inv_node( fo );
            }
          } );
          DEBUG_PRINT_MIG( ntk );
          minimized = true;
        }
      } );
      break;
    }
  }

private:
  int two_level( node<Ntk> n )
  {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) || ntk.is_dead( n ) )
      return 0;
    int gain = one_level( n );
    ntk.foreach_fanout( n, [&]( auto const& fo ) {
      int gain_fo = one_level( fo );
      if ( is_fanout_comp( n, fo ) )
        gain_fo -= 2;
      else
        gain_fo += 2;
      INFO( "          gain_fo: " << gain_fo );
      if ( gain_fo > 0 )
        gain += gain_fo;
    } );
    INFO( "        2 level gain: " << gain );
    return gain;
  }

  int one_level( node<Ntk> n )
  {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) || ntk.is_dead( n ) )
      return 0;

    // number of complemented inputs
    uint32_t num_compl_inp = 0;
    uint32_t num_non_compl_inp = 0;
    ntk.foreach_fanin( n, [&]( auto const& fi ) {
      if ( ntk.is_constant( fi.index ) )
        return;
      if ( fi.complement )
        num_compl_inp++;
      else
        num_non_compl_inp++;
    } );
    DEBUG( "    num_compl_inp: " << num_compl_inp );
    DEBUG( "    num_non_compl_inp: " << num_non_compl_inp );

    // number of complemented outputs
    uint32_t num_compl_out = 0;
    uint32_t num_non_compl_out = 0;
    ntk.foreach_fanout( n, [&]( auto const& fo ) {
      if ( is_fanout_comp( n, fo ) )
        num_compl_out++;
      else
        num_non_compl_out++;
    } );
    DEBUG( "    num_compl_out: " << num_compl_out );
    DEBUG( "    num_non_compl_out: " << num_non_compl_out );

    // number of complemented po
    uint32_t num_compl_po = 0;
    uint32_t num_non_compl_po = 0;
    ntk.foreach_po( [&num_compl_po, &num_non_compl_po, &n]( auto const& po ) {
      if ( po.index == n )
      {
        if ( po.weight )
          num_compl_po++;
        else
          num_non_compl_po++;
      }
    } );

    // total gain
    int gain = num_compl_inp + num_compl_out + num_compl_po - num_non_compl_inp - num_non_compl_out - num_non_compl_po;
    INFO( "        gain: " << gain );
    return gain;
  }

  node<Ntk> inv_node( node<Ntk> n, bool have_comp_out = true )
  {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return n;

    auto a = ntk._storage->nodes[n].children[0];
    auto b = ntk._storage->nodes[n].children[1];
    auto c = ntk._storage->nodes[n].children[2];

    a.weight ^= 1;
    b.weight ^= 1;
    c.weight ^= 1;

    auto new_node = !create_maj_without_changing_comp( a, b, c );
    ntk.replace_in_outputs( n, new_node );
    ntk.foreach_fanout( n, [&]( auto const& fo ) {
      if ( have_comp_out || is_fanout_comp( n, fo ) )
        ntk.replace_in_node( fo, n, new_node );
    } );
    if ( ntk.fanout_size( n ) == 0 )
      ntk.take_out_node( n );
    return ntk.get_node( new_node );
  }

  bool is_fanout_comp( node<Ntk> in, node<Ntk> out )
  {
    bool ret = false;
    ntk.foreach_fanin( out, [&in, &ret]( auto const& fi ) {
      if ( fi.index == in )
        ret = fi.complement != 0;
    } );
    return ret;
  }

  signal<Ntk> create_maj_without_changing_comp( signal<Ntk> a, signal<Ntk> b, signal<Ntk> c )
  {
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }
    else
    {
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }

    /* trivial cases */
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? a : c;
    }
    else if ( b.index == c.index )
    {
      return ( b.complement == c.complement ) ? b : a;
    }

    typename Ntk::storage::element_type::node_type new_node;
    new_node.children[0] = a;
    new_node.children[1] = b;
    new_node.children[2] = c;

    /* structural hashing */
    const auto it = ntk._storage->hash.find( new_node );
    if ( it != ntk._storage->hash.end() )
    {
      return { it->second, 0 };
    }

    const auto index = ntk._storage->nodes.size();

    if ( index >= .9 * ntk._storage->nodes.capacity() )
    {
      ntk._storage->nodes.reserve( static_cast<uint64_t>( 3.1415f * index ) );
      ntk._storage->hash.reserve( static_cast<uint64_t>( 3.1415f * index ) );
    }

    ntk._storage->nodes.push_back( new_node );

    ntk._storage->hash[new_node] = index;

    /* increase ref-count to children */
    ntk._storage->nodes[a.index].data[0].h1++;
    ntk._storage->nodes[b.index].data[0].h1++;
    ntk._storage->nodes[c.index].data[0].h1++;

    for ( auto const& fn : ntk._events->on_add )
    {
      ( *fn )( index );
    }

    return { index, 0 };
  }

private:
  Ntk& ntk;
  mig_inv_minimization_params const& ps;
  mig_inv_minimization_stats& st;
};
} // namespace detail

/*! \brief Majority algebraic depth rewriting.
 *
 * This algorithm tries to rewrite a network with majority gates for depth
 * optimization using the associativity and distributivity rule in
 * majority-of-3 logic.  It can be applied to networks other than MIGs, but
 * only considers pairs of nodes which both implement the majority-of-3
 * function.
 *
 * **Required network functions:**
 * - `get_node`
 * - `level`
 * - `update_levels`
 * - `create_maj`
 * - `substitute_node`
 * - `foreach_node`
 * - `foreach_po`
 * - `foreach_fanin`
 * - `is_maj`
 * - `clear_values`
 * - `set_value`
 * - `value`
 * - `fanout_size`
 *
   \verbatim embed:rst

  .. note::

      The implementation of this algorithm was heavily inspired by an
      implementation from Luca Amarù.
   \endverbatim
 */
template<class Ntk>
void mig_inv_minimization( Ntk& ntk, mig_inv_minimization_params const& ps = {}, mig_inv_minimization_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_replace_in_outputs_v<Ntk>, "Ntk does not implement the replace_in_outputs method" );
  static_assert( has_replace_in_node_v<Ntk>, "Ntk does not implement the replace_in_node method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

  mig_inv_minimization_stats st;
  fanout_view<Ntk> fo_ntk( ntk );
  detail::mig_inv_minimization_impl<fanout_view<Ntk>> p( fo_ntk, ps, st );
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle