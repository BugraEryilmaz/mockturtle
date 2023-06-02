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
  \file mig_inv_propogation.hpp
  \brief MIG Fanout Restriction Algorithm

  \author Bugra Eryilmaz
*/

#pragma once

#include "../networks/mig.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/fanout_view.hpp"

#include <iostream>
#include <optional>
#include <queue>

#define DLEVEL 2
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

/*! \brief Parameters for mig_inv_propogation.
 *
 * The data structure `mig_inv_propogation_params` holds configurable
 * parameters with default arguments for `mig_inv_propogation`.
 */
struct mig_inv_propogation_params
{
};

/*! \brief Statistics for mig_inv_propogation. */
struct mig_inv_propogation_stats
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
class mig_inv_propogation_impl
{
public:
  mig_inv_propogation_impl( Ntk& ntk, mig_inv_propogation_params const& ps, mig_inv_propogation_stats& st )
      : ntk( ntk ), ps( ps ), st( st )
  {
  }
  void run()
  {
    ntk.foreach_po( [&]( auto const& f ) {
      if ( f.weight )
      {
        q.push( f.index );
      }
      else
      {
        q.push( f.index );
      }
    } );

    while ( !q.empty() )
    {
      auto n = q.front();
      q.pop();
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_dead( n ) )
        continue;
      if ( has_comp_fanout( n ) )
      {
        INFO( "Visiting node " << n << " has complement fanout" );
        inv_node( n );
        DEBUG_PRINT_MIG( ntk );
      }
      else
      {
        INFO( "Visiting node " << n << " has no complement fanout" );
      }
      ntk.foreach_fanin( n, [&]( auto const& fi ) {
        q.push( fi.index );
      } );
    }
  }

private:
  bool has_comp_fanout( node<Ntk> n )
  {
    bool ret = false;
    ntk.foreach_gate( [&]( auto const& fo ) {
      DEBUG( "Fanout of " << n << " is " << fo << " is comp: " << is_fanout_comp( n, fo ) );
      if ( is_fanout_comp( n, fo ) )
        ret = true;
    } );
    if ( ret )
      return true;
    ntk.foreach_po( [&]( auto const& po ) {
      if ( po.index == n && po.weight )
        ret = true;
    } );
    return ret;
  }
  node<Ntk> inv_node( node<Ntk> n, bool have_comp_out = false )
  {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return n;

    auto a = ntk._storage->nodes[n].children[0];
    auto b = ntk._storage->nodes[n].children[1];
    auto c = ntk._storage->nodes[n].children[2];

    a.weight ^= 1;
    b.weight ^= 1;
    c.weight ^= 1;

    print_all_fanout();
    auto new_node = !create_maj_without_changing_comp( a, b, c );
    print_all_fanout();
    replace_in_outputs_cond_comp( n, new_node, have_comp_out );
    ntk.foreach_gate( [&]( auto const& fo ) {
      if ( have_comp_out || is_fanout_comp( n, fo ) )
      {
        auto ret = ntk.replace_in_node( fo, n, new_node );
        if ( ret )
        {
          char const* const comp = ret->second.complement ? "!" : "";
          INFO( "Replaced node " << fo << " with " << comp << ret->second.index );
          ntk.substitute_node( fo, ret->second );
        }
      }
    } );
    if ( ntk.fanout_size( n ) == 0 )
      ntk.take_out_node( n );
    return ntk.get_node( new_node );
  }
  void print_all_fanout()
  {
    for ( int i = 0; i < ntk._fanout.size(); i++ )
    {
      auto f = ntk._fanout[i];
      std::cout << "Fanout of " << i << " is: ";
      for ( auto fo : f )
      {
        std::cout << fo << " ";
      }
      std::cout << std::endl;
    }
  }

  void replace_in_outputs_cond_comp( node<Ntk> const& old_node, signal<Ntk> const& new_signal, bool have_comp_out = false )
  {
    if ( ntk.is_dead( old_node ) )
      return;

    for ( auto& output : ntk._storage->outputs )
    {
      if ( output.index == old_node && ( have_comp_out || output.weight ) )
      {
        output.index = new_signal.index;
        output.weight ^= new_signal.complement;

        if ( old_node != new_signal.index )
        {
          // increment fan-out of new node
          ntk._storage->nodes[new_signal.index].data[0].h1++;
          // decrement fan-out of old node
          ntk._storage->nodes[old_node].data[0].h1--;
        }
      }
    }
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
  mig_inv_propogation_params const& ps;
  mig_inv_propogation_stats& st;
  std::queue<node<Ntk>> q;
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
void mig_inv_propogation( Ntk& ntk, mig_inv_propogation_params const& ps = {}, mig_inv_propogation_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_replace_in_outputs_v<Ntk>, "Ntk does not implement the replace_in_outputs method" );
  static_assert( has_replace_in_node_v<Ntk>, "Ntk does not implement the replace_in_node method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

  mig_inv_propogation_stats st;
  fanout_view<Ntk> fo_ntk( ntk );
  detail::mig_inv_propogation_impl<fanout_view<Ntk>> p( fo_ntk, ps, st );
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} // namespace mockturtle