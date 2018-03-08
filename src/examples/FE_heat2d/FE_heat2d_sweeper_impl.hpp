
//#include "FE_sweeper.hpp"

#include "assemble2d.hpp"

/*use grid-glue */
#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>

#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/gridglue.hh>
/*end grid-glue */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include <vector>
using std::shared_ptr;
using std::vector;

//#include <leathers/push>
//#include <leathers/all>
#include <boost/math/constants/constants.hpp>
//#include <leathers/pop>
using boost::math::constants::pi;
using boost::math::constants::two_pi;
using boost::math::constants::pi_sqr;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>

namespace Helper
{
  template <class GridView>
  struct VerticalFaceDescriptor
  : public Dune::GridGlue::ExtractorPredicate<GridView,1>
  {
    virtual bool contains(const typename GridView::template Codim<0>::Entity& element,
			  unsigned int face) const
    {
      const int dim = GridView::dimension;
      const auto& refElement = Dune::ReferenceElements<double, dim>::general(element.type());
      
      // Number of corners of the element face
      int numVertices = refElement.size(face, 1, dim);
      
      for (int i=0; i<numVertices; i++)
	if ( std::abs(element.geometry().corner(refElement.subEntity(face,1,i,dim))[0] ) > 1e-6 )
	  return false;
	
	return true;
    }
  };
}

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      template<class SweeperTrait, typename Enabled>
      void
      Heat2d_FE<SweeperTrait, Enabled>::init_opts()
      {
#if 0
        config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");
#endif
      }


        template<class SweeperTrait, typename Enabled>
      Heat2d_FE<SweeperTrait, Enabled>::Heat2d_FE(const size_t nelements, const size_t basisorder, const int finer)
        :   IMEX<SweeperTrait, Enabled>(),
        lower_mv({-0.5, 2.0}), upper_mv({0.0, 2.5}),
        hgtmv(lower_mv, upper_mv, std::array<int, 2>({2, 2})),
        grid_mv(new GridType_MV(hgtmv, mf)),
        v0(5.0), alpha(1.0e-3), sourceVal(0.0), vel(-0.1) //vel(-0.1)
      {

        //const auto dim = 2;//SweeperTrait::DIM;        

        this->grid        = createGrid(nelements);


        grid->globalRefine(finer);
        this->finer = finer;

        typedef GridType::LeafGridView GridView;
        //typedef GridType::LevelGridView LevelView;

        GridType::LeafGridView gridView       = grid->leafGridView();
        //GridType::LevelGridView gridView_coase = grid_coarse->levelGridView(1);

        std::cout << "***** Anzahl der finiten Elemente " << nelements << std::endl;
        //std::cout << "***** Ordnung der Basis " << SweeperTrait::BASE_ORDER << std::endl;

        this->basis       = std::make_shared<BasisFunction>(gridView);
        //this->basis_coase = std::make_shared<BasisFunction>(gridView_coase);

        const auto bs = basis->size();
        std::cout << "***** Basis erstellt mit " <<  basis->size() << " Elementen " << std::endl;

        this->encap_factory()->set_size(bs);




      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename Heat2d_FE<SweeperTrait, Enabled>::GridType> 
      Heat2d_FE<SweeperTrait, Enabled>::createGrid(unsigned nelements)
      {
	 Dune::FieldVector<double,dim> h = {1,4};
        array<int,dim> n;
	n[0] = nelements;
	n[1] = 4*nelements; 
	std::shared_ptr<GridType > ret (make_shared<GridType>(h,n));
	return ret;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat2d_FE<SweeperTrait, Enabled>::set_options()
      {


        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        std::cout << "nu    "   << this->_nu << std::endl;

        int num_nodes = this->get_quadrature()->get_num_nodes();

        assembleProblem(basis, this->A_dune, this->M_dune/*, b_dune*/);


      }

        template<class SweeperTrait, typename Enabled>
        template<typename Basis>
      void
      Heat2d_FE<SweeperTrait, Enabled>::assemble(Basis &basis){

        assembleProblem(basis, this->A_dune, this->M_dune /*, b_dune*/);


      };

      template<class SweeperTrait, typename Enabled>
      void
      Heat2d_FE<SweeperTrait, Enabled>::assemble(){

        assembleProblem(basis, this->A_dune, this->M_dune /*, b_dune*/);


      };




      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat2d_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();
        result->data()=0.0;
        #if 0
        
        //const auto dim = 2;//SweeperTrait::DIM;
        spatial_t nu = this-> _nu;
        
        auto exact_solution = [t, nu, dim](const FieldVector<double,dim>&x){
        double solution=1.0;
        for(int i=0; i<dim; i++){solution *= std::sin(pi<spatial_t>() * x[i]);}
        return solution * std::exp(-t * dim * pi_sqr<spatial_t>() * nu);
      };
      
      auto N_x = [t](const FieldVector<double,dim>&x){
      return x;
      
      };
      
      BlockVector<FieldVector<double,dim>> x_node;
      interpolate(*basis, x_node, N_x);
      
      interpolate(*basis, result->data(), exact_solution);
      #endif
      
      return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat2d_FE<SweeperTrait, Enabled>::initial()
      {
        auto result = this->get_encap_factory().create();
        result->data()=1.0;
#if 0
        for(int i(0); i < result->data().size(); ++i)
          result->data()[i] = i;
#endif
        return result;
      }
      
      template<class SweeperTrait, typename Enabled>
      void
      Heat2d_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat2d_FE<SweeperTrait, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" << (t + dt * nodes[m])
                     << "      |abs residual| = " << this->_abs_res_norms[m]
                     << "      |rel residual| = " << this->_rel_res_norms[m]
//                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
//                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" << (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " << this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " << this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat2d_FE<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat2d_FE<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }


      //typedef Dune::YaspGrid<2> GridType;
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename Heat2d_FE<SweeperTrait, Enabled>::GridType >
      Heat2d_FE<SweeperTrait, Enabled>::get_grid() const
      {
        return grid;
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat2d_FE<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat2d_FE<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                            const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

#if 0
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat2d_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        //ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating EXPLICIT part at t=" << t);	

        auto result = this->get_encap_factory().create();	
        result->zero();

	auto swap = this->get_encap_factory().create();
	
	//swap->data() = b_dune;
	swap->zero();
	//fastGrid(t, u->data(), swap->data());
	fastRect(t, u->data(), swap->data());
	MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dune);

        SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dune,1.0);

        CGSolver<VectorType> cg(linearOperator,
                                preconditioner,
                                1e-10, // desired residual reduction factor
                                500,    // maximum number of iterations
                                0);    // verbosity of the solver
        InverseOperatorResult statistics ;

        cg.apply(result->data(), swap->data() , statistics );

        this->_num_expl_f_evals++;

        return result;
      }
#endif

      template<class SweeperTrait, typename Enabled>
      void Heat2d_FE<SweeperTrait, Enabled>::fastGrid(double t, const VectorType& yIn, VectorType& out) const
      {	
        using namespace Dune;
	//const double center = 2.25-0.1*t;
	const double dY = vel*t; //-0.1*t;
	//move grid_mv by dY (indirectly through geometrygrid)
	mf.dx[1] = dY;
	
	typedef typename GridType_MV::LeafGridView GridView_MV;
	
	Helper::VerticalFaceDescriptor<GridView> facePredicate0;
	Helper::VerticalFaceDescriptor<GridView_MV> facePredicate1;
	
	typedef Dune::GridGlue::Codim1Extractor<GridView> Extractor0;
	typedef Dune::GridGlue::Codim1Extractor<GridView_MV> Extractor1;
	
	Extractor0 domEx(grid->leafGridView(), facePredicate0);
	Extractor1 tarEx(grid_mv->leafGridView(), facePredicate1);
	
	typedef Dune::GridGlue::GridGlue<Extractor0,Extractor1> GlueType;
	
	// Backend for the computation of the remote intersections
	Dune::GridGlue::ContactMerge<dim,double> merger;
	GlueType glue(domEx, tarEx, &merger);
	
	glue.build();
	
	const GridView::IndexSet& indexSet0 = grid->leafGridView().indexSet();
	const typename GridView_MV::IndexSet& indexSet1 = grid_mv->leafGridView().indexSet();
	
	typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> TestFECache;
	typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> FiniteElementCache0;
	FiniteElementCache0 cache0, cache1;
	TestFECache testCache;
	
	for (const auto& intersection : intersections(glue))
	{
	  const FiniteElementCache0::FiniteElementType& nonmortarFiniteElement = cache0.get(intersection.inside().type());
	  const FiniteElementCache0::FiniteElementType& mortarFiniteElement    = cache1.get(intersection.outside().type());
	  const TestFECache::FiniteElementType&         testFiniteElement      = testCache.get(intersection.inside().type());
	  
	  // Select a quadrature rule:  Use order = 2 just for simplicity
	  int quadOrder = 2;
	  const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), quadOrder);
	  
	  // Loop over all quadrature points
	  for (size_t l=0; l<quad.size(); l++)
	  {
	    // compute integration element of overlap
	    double integrationElement = intersection.geometry().integrationElement(quad[l].position());
	    
	    // quadrature point positions on the reference element
	    FieldVector<double,dim> nonmortarQuadPos = intersection.geometryInInside().global(quad[l].position());
	    //FieldVector<double,dim> mortarQuadPos    = intersection.geometryInOutside().global(quad[l].position());
	    
	    //evaluate all shapefunctions at the quadrature point
	    std::vector<FieldVector<double,1> > nonmortarValues,testValues; // mortarValues
	    
	    nonmortarFiniteElement.localBasis().evaluateFunction(nonmortarQuadPos,nonmortarValues);
	    //mortarFiniteElement   .localBasis().evaluateFunction(mortarQuadPos,mortarValues);
	    testFiniteElement     .localBasis().evaluateFunction(nonmortarQuadPos,testValues);
	    
	    double uh(0.0);
	    for (size_t j=0; j<nonmortarValues.size(); j++) { 
	      auto r = indexSet0.subIndex(intersection.inside(), j, dim);
	      uh += nonmortarValues[j]*yIn[r];
	    }
	    const double fVal = (v0-uh)*alpha+sourceVal;
	    // Loop over all shape functions of the test space
	    for (size_t j=0; j<testFiniteElement.size(); j++)
	    {
	      int testIdx = indexSet0.subIndex(intersection.inside(),j,dim);
	      out[testIdx] += integrationElement*quad[l].weight()*testValues[j]*fVal;
	      
	    }
	    
	  }
	  
	}
      }
      
      template<class SweeperTrait, typename Enabled>
      void Heat2d_FE<SweeperTrait, Enabled>::fastRect(double t, const VectorType& yIn, VectorType& out) const
      {
        const double center = 2.25-0.1*t;
        auto rect = [center](const auto& x) { double d=(center-x[1])/0.25; return  std::abs(d) <= 1.0 ? cos(d*M_PI/2) : 0.0 ;  };
        //auto rect = [center](const auto& x) { double d=(center-x[1])/0.25; return  std::abs(d) <= 1.0 ? 1.0 : 0.0 ;  };
        //add the neumann part
        fastBoundary(yIn, [&](double uh, const auto& posGlobal) { return rect(posGlobal)*alpha*(v0 /*-uh*/); }, out);
      }
      
      template<class SweeperTrait, typename Enabled>
      template<typename Flux >
      void Heat2d_FE<SweeperTrait, Enabled>::fastBoundary(const VectorType& yIn, const Flux& flux, VectorType& out) const
      {
        using namespace Dune;

        auto localView = basis->localView();
        auto localIndexSet = basis->localIndexSet();
        auto gridView = grid->leafGridView();
        for (const auto& bdEl : elements(gridView))
        {
          localView.bind(bdEl);
          localIndexSet.bind(localView);
          for(const auto& inter : intersections(gridView, bdEl))
          {
            //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
            if(inter.boundary() && (std::abs(inter.geometry().center()[0]) < 1.0e-10)) {
              //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
              const int qOrder = 2*dim+1+3; // 2;
              auto quadRule = QuadratureRules<double, dim-1>::rule(inter.type(), qOrder);
              const auto& localFiniteElement = localView.tree().finiteElement();
              for(const auto& qPos : quadRule) {
                //std::cout << "qPos:" << qPos.position(); // <<std::endl;
                const double det = inter.geometry().integrationElement(qPos.position());
                std::vector<FieldVector<double,1> > shapeFunctionValues;
                //localFiniteElement.localBasis().evaluateFunction(inter.geometryInInside().global(qPos.position()), shapeFunctionValues);
                //std::cout << "\ninInside.type():" <<inter.geometryInInside().type() ;
                //std::cout << "\nqGlobal:" << inter.geometryInInside().global(qPos.position()) ;
                auto posInRef = inter.geometryInInside().global(qPos.position());
                localFiniteElement.localBasis().evaluateFunction(posInRef, shapeFunctionValues); //ist das richtig?? ist das global in Ref, also aus sicht der Linie global
                auto posGlobal = bdEl.geometry().global(posInRef);
                //std::cout << "globalPosBoundary:" << posGlobal << std::endl;
                double uh = 0.0;
                for(unsigned i(0); i < localFiniteElement.size(); ++i) {
                  const auto row = localIndexSet.index(i);
                  uh += shapeFunctionValues[i]*yIn[row];
                }
                const double fVal = flux(uh, posGlobal); 
                //std::cout << "shapeVals:[";
                for(unsigned i(0); i < localFiniteElement.size(); ++i) {
                  //std::cout << shapeFunctionValues[i] << " ";
                  const auto row = localIndexSet.index(i);
                  out[row] += qPos.weight()*fVal*shapeFunctionValues[i]*det;
                }
                //std::cout << std::endl;
              }
              
            }
          }
        }
      }
      
      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat2d_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {

        using namespace Dune;

        //ML_CVLOG(4, this->get_logger_id(), LOG_FIXED << "evaluating IMPLICIT part at t=" << t);


        auto result = this->get_encap_factory().create();
        //auto rhs = this->get_encap_factory().create();
        //auto rhs2 = this->get_encap_factory().create();
        double nu =this->_nu;
        //u->data()*= nu;
        result->data().resize(u->data().size());
        this->A_dune.mv(u->get_data(), result->data());	
        //A_dune.mmv(u->get_data(), rhs2->data());



        //auto DirichletValues = [] (auto x) {return (x[0]<1e-8 or x[0]>0.9999) ? 0 : x[0];};
#if 0

        Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dune,1.0);

        CGSolver<VectorType> cg(linearOperator,
                                preconditioner,
                                1e-10, // desired residual reduction factor
                                500,    // maximum number of iterations
                                0);    // verbosity of the solver
        InverseOperatorResult statistics ;

        cg.apply(result->data(), rhs->data() , statistics );

	/*
        auto var = this->get_encap_factory().create();
        M_dune.mv(result->data(), var->data());
        auto neu = this->get_encap_factory().create();
        rhs2->scaled_add(-1.0 , var)  ;
	*/

#endif
        result->data() *= -nu;	

        //add the source
        fastRect(t, u->get_data(), result->data());
        return result;


      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat2d_FE<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {
        
        using namespace Dune;
        
        Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune, swap ;
        M_rhs_dune.resize(rhs->get_data().size());
        swap.resize(rhs->get_data().size());
        
        M_rhs_dune = rhs->data();
        //std::cout << "norm(rhs): " << rhs->norm0() << std::endl;
        //M_dune.mv(rhs->data(), M_rhs_dune);
        //swap = M_rhs_dune;
        
        //add the source
        swap=0.0; fastRect(t, u->get_data(), swap);
        //M_rhs_dune += dt*swap;
        M_rhs_dune.axpy(dt, swap);
        
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->A_dune);
        M_dtA_dune *= (dt * this->_nu);
        M_dtA_dune += this->M_dune;
        
        MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);
        
        SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);
        
        CGSolver<VectorType> cg(linearOperator,
                                preconditioner,
                                1e-10, // desired residual reduction factor
                                500,    // maximum number of iterations
                                10);    // verbosity of the solver
        
        InverseOperatorResult statistics ;
        
        cg.apply(u->data(), M_rhs_dune , statistics ); //rhs ist nicht constant!!!!!!!!!
        
        
        #if 0
        ML_CVLOG(4, this->get_logger_id(),
        LOG_FIXED << "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);
        #endif
        #if 0
        swap = u->data();
        swap -= rhs->data();
        swap *= (1.0/dt);
        
        MatrixAdapter<MatrixType,VectorType,VectorType> linearOperatorM(M_dune);
        
        SeqILU0<MatrixType,VectorType,VectorType> preconditionerM(M_dune,1.0);
        
        CGSolver<VectorType> cgM(linearOperatorM,
        preconditionerM,
        1e-10, // desired residual reduction factor
        500,    // maximum number of iterations
        0);    // verbosity of the solver
        
        InverseOperatorResult statisticsM ;
        
        cgM.apply(f->data(), swap, statisticsM);
        #else
        /*
        M_dune.mv(u->get_data(), f->data());
        f->data()-=rhs->data();
        f->data() *= (1.0/dt);
        */
        this->A_dune.mv(u->data(), f->data()); f->data()*=-this->_nu;
        f->data().axpy(1.0, swap);
        //std::cout << "u after solve: " ; for(auto& d : u->data()) std::cout << " " << d; std::cout << std::endl;
        
        /*for (size_t i = 0; i < u->data().size(); i++) {
         *   f->data()[i] = 0.0;(u->data()[i] - rhs->data()[i]) / (dt);
      }*/
        #endif
        
        
        this->_num_impl_solves++;
        
        
        
        
        
      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
}  // ::pfasst
