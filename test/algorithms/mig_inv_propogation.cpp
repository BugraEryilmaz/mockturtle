#include <catch.hpp>

#include <mockturtle/algorithms/mig_inv_propogation.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/views/fanout_view.hpp>

using namespace mockturtle;

namespace mockturtle::mig_prop
{
uint64_t complement_count( mig_network const& mig )
{
  uint64_t count{ 0 };
  mig.foreach_node( [&]( auto n ) {
    mig.foreach_fanin( n, [&]( auto const& f ) {
      if ( f.complement )
      {
        count++;
      }
    } );
  } );
  mig.foreach_po( [&count]( auto const& po ) {
    if ( po.weight )
      count++;
  } );

  return count;
}

uint64_t complement_count_except_pi( mig_network const& mig )
{
  uint64_t count{ 0 };
  mig.foreach_node( [&]( auto n ) {
    mig.foreach_fanin( n, [&count, &mig]( auto const& f ) {
      if ( f.complement && !mig.is_pi( f.index ) && !mig.is_constant( f.index ) )
      {
        count++;
      }
    } );
  } );
  mig.foreach_po( [&count, &mig]( auto const& po ) {
    if ( po.weight && !mig.is_pi( po.index ) && !mig.is_constant( po.index ) )
      count++;
  } );

  return count;
}
} // namespace mockturtle::mig_prop
TEST_CASE( "MIG inverter propogation one level", "[mig_inv_propogation]" )
{
  mig_network mig;

  const auto a = mig.create_pi();
  const auto b = mig.create_pi();
  const auto c = mig.create_pi();
  const auto d = mig.create_pi();
  const auto e = mig.create_pi();

  const auto f1 = mig.create_maj( a, b, !c );                       // 1
  const auto f2 = mig.create_maj( a, b, mig.get_constant( true ) ); // 0
  const auto f3 = mig.create_maj( !f1, f2, a );
  const auto f4 = mig.create_maj( !f1, f2, b );
  const auto f5 = mig.create_maj( !f1, f2, c );
  const auto f6 = mig.create_maj( f1, !f2, a ); // 1
  const auto f7 = mig.create_maj( f4, !f2, b );
  const auto f8 = mig.create_maj( f5, !f2, c );
  const auto f9 = mig.create_maj( f6, !f2, d );
  const auto f10 = mig.create_maj( f7, !f2, e );

  mig.create_po( f3 );
  mig.create_po( f4 );
  mig.create_po( f5 );
  mig.create_po( !f6 );
  mig.create_po( f7 );
  mig.create_po( f8 );
  mig.create_po( !f9 );
  mig.create_po( !f10 );

  print( mig );

  mig_inv_propogation_params pr;
  mig_inv_propogation_stats st;
  mig_inv_propogation( mig, pr, &st );

  print( mig );

  CHECK( mig_prop::complement_count_except_pi( mig ) == 0 );
}

// TEST_CASE( "MIG inverter propogation two level", "[mig_inv_propogation]" )
// {
//   mig_network mig;

//   const auto a = mig.create_pi();
//   const auto b = mig.create_pi();
//   const auto c = mig.create_pi();

//   const auto f1 = mig.create_maj( a, b, !c );
//   const auto f2 = mig.create_maj( a, b, c );
//   const auto f3 = mig.create_maj( !f1, f2, a );
//   const auto f4 = mig.create_maj( !f3, b, c );
//   const auto f5 = mig.create_maj( f1, b, !a );
//   const auto f6 = mig.create_maj( !f1, f2, b );

//   mig.create_po( f4 );
//   mig.create_po( !f5 );
//   mig.create_po( f6 );

//   print( mig );

//   mig_inv_propogation_params pr;
//   mig_inv_propogation_stats st;
//   mig_inv_propogation( mig, pr, &st );

//   print( mig );

//   CHECK( mig_prop::complement_count_except_pi( mig ) == 0 );
// }
