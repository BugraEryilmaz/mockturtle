#include <catch.hpp>

#include <mockturtle/algorithms/mig_inv_minimization.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/debugging_utils.hpp>
#include <mockturtle/views/fanout_view.hpp>

using namespace mockturtle;

namespace mockturtle::mig_comp
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
} // namespace mockturtle::mig_comp
TEST_CASE( "MIG inverter minimization one level", "[mig_inv_minimization]" )
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
  mig.create_po( f9 );
  mig.create_po( f10 );

  print( mig );

  uint32_t gates_old = mig.num_gates();
  uint64_t comp_count_old = mig_comp::complement_count( mig );

  mig_inv_minimization_params pr;
  mig_inv_minimization_stats st;
  mig_inv_minimization( mig, pr, &st );

  print( mig );

  CHECK( mig.num_gates() == gates_old );

  uint64_t comp_count = mig_comp::complement_count( mig );

  CHECK( comp_count_old > comp_count );
  CHECK( comp_count_old - comp_count == st.num_inverters_removed );
  CHECK( st.num_inverters_removed == 2 );
}

TEST_CASE( "MIG inverter minimization two level", "[mig_inv_minimization]" )
{
  mig_network mig;

  const auto a = mig.create_pi();
  const auto b = mig.create_pi();
  const auto c = mig.create_pi();
  const auto d = mig.create_pi();
  const auto e = mig.create_pi();

  const auto f1 = mig.create_maj( a, b, !c );
  const auto f2 = mig.create_maj( a, b, c );
  const auto f3 = mig.create_maj( !f1, f2, a );
  const auto f4 = mig.create_maj( !f3, b, c );
  const auto f5 = mig.create_maj( f1, b, !a );
  const auto f6 = mig.create_maj( !f1, f2, b );

  mig.create_po( f4 );
  mig.create_po( !f5 );
  mig.create_po( f6 );

  print( mig );

  uint32_t gates_old = mig.num_gates();
  uint64_t comp_count_old = mig_comp::complement_count( mig );

  mig_inv_minimization_params pr;
  mig_inv_minimization_stats st;
  mig_inv_minimization( mig, pr, &st );

  print( mig );

  CHECK( mig.num_gates() == gates_old );

  uint64_t comp_count = mig_comp::complement_count( mig );

  CHECK( comp_count_old > comp_count );
  CHECK( comp_count_old - comp_count == st.num_inverters_removed );
  CHECK( st.num_inverters_removed == 2 );
}
