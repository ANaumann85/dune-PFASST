//#include <dune/istl/matrixmarket.hh>
//#include <leathers/push>
//#include <leathers/all>

#include <dune/common/function.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include <dune/istl/matrixmarket.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <vector>

//#include <leathers/pop>




template<class LocalView, class MatrixType>
void assembleElementA(const LocalView &localView, MatrixType &elementMatrix) {

  using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();

  elementMatrix.setSize(localFiniteElement.size(), localFiniteElement.size());
  elementMatrix = 0;




  int order = 2 * (dim * localFiniteElement.localBasis().order() - 1);


  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);


  for (size_t pt = 0; pt < quad.size(); pt++) {   // Loop over all quadrature points


    const auto quadPos = quad[pt].position(); // Position of the current quadrature point in the reference element

    const auto jacobian = geometry.jacobianInverseTransposed(quadPos); // The transposed inverse Jacobian of the map from the reference element to the element

    const auto integrationElement = geometry.integrationElement(quadPos); // The multiplicative factor in the integral transformation formula


    std::vector<FieldMatrix<double, 1, dim> > referenceGradients;
    localFiniteElement.localBasis().evaluateJacobian(quadPos, referenceGradients); // The gradients of the shape functions on the reference element

    std::vector<FieldVector<double, dim> > gradients(referenceGradients.size());
    for (size_t i = 0; i < gradients.size(); i++)
      jacobian.mv(referenceGradients[i][0], gradients[i]);     // Compute the shape function gradients on the real element


    for (size_t i = 0; i < elementMatrix.N(); i++)
      for (size_t j = 0; j < elementMatrix.M(); j++)
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                += (gradients[i] * gradients[j]) * quad[pt].weight() * integrationElement;

  }
}

template<class LocalView, class MatrixType>
void assembleElementM(const LocalView &localView, MatrixType &elementMatrix) {

  using namespace Dune;
  using Element = typename LocalView::Element;

  auto element = localView.element();
  const int dim = Element::dimension;
  auto geometry = element.geometry();

  const auto &localFiniteElement = localView.tree().finiteElement();


  elementMatrix.setSize(localFiniteElement.size(), localFiniteElement.size());
  elementMatrix = 0;


  int order = 2 * (dim * localFiniteElement.localBasis().order() );

  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  for (size_t pt = 0; pt < quad.size(); pt++) {

    const auto quadPos = quad[pt].position();

    //const auto transformationsfunktion = geometry.

    const auto integrationElement = geometry.integrationElement(quadPos);


    std::vector<FieldVector<double, 1> > shapeFunctionValues;
    localFiniteElement.localBasis().evaluateFunction(quadPos , shapeFunctionValues);//element.geometryInInside().global(quadPos);

    /*for (size_t i = 0; i < shapeFunctionValues.size(); i++) {
      std::cout << shapeFunctionValues[i] << std::endl;
    }*/
    //shapeFunctionValues[i] = element.geometry().global(shapeFunctionValues[i]);}


    for (size_t i = 0; i < elementMatrix.N(); i++){
      for (size_t j = 0; j < elementMatrix.M(); j++) {
        elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
                += (shapeFunctionValues[i] * shapeFunctionValues[j]) * quad[pt].weight() * integrationElement;

        //std::cout << "der eintrag der element matrix " << elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)] << std::endl;
      }
    }
    /*for (size_t i = 0; i < elementMatrix.N(); i++){
      for (size_t j = 0; j < elementMatrix.M(); j++) {

        std::cout << "der eintrag der element matrix " << elementMatrix[localView.tree().localIndex(i)][localView.tree().localIndex(j)] << std::endl;
        std::cout << "quad[pt].weight() " <<  quad[pt].weight() << std::endl;
        std::cout << "integrationElement " <<  integrationElement << std::endl;
        std::cout << "integrationElement " <<  integrationElement << std::endl;
      }
    }*/

    //std::exit(0);

  }
}

template<typename Flux, class Basis, class GridView, typename VT>
void fastBoundary(Flux& flux, const Basis& basis, GridView& gridView, VT& out)
{
  using namespace Dune;
  static const unsigned dim = 2;
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();
  for (const auto& bdEl : elements(gridView))
  {
    localView.bind(bdEl);
    localIndexSet.bind(localView);
    for(const auto& inter : intersections(gridView, bdEl))
    {
      //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
      if(inter.boundary() && (std::abs(inter.geometry().center()[0]) < 1.0e-10)) {
        //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
        const int qOrder = 2;
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
//#pragma message("###assemble2d.hpp:174: constant source###")
            //uh += shapeFunctionValues[i]*yIn[row];
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


template<class Basis>
void getOccupationPattern(const Basis &basis, Dune::MatrixIndexSet &nb) {
  using namespace Dune;
  nb.resize(basis->size(), basis->size());
  auto gridView = basis->gridView();
  auto localView = basis->localView();
  auto localIndexSet = basis->localIndexSet();
  for (const auto &element : elements(gridView)) { // A loop over all elements of the grid

    localView.bind(element);
    localIndexSet.bind(localView);

    for (size_t i = 0; i < localIndexSet.size(); i++) {


      auto row = localIndexSet.index(i);
      for (size_t j = 0; j < localIndexSet.size(); j++) {

        auto col = localIndexSet.index(j);
        nb.add(row, col);
      }
    }
  }

}



template<class Basis>
void assembleProblem(const Basis &basis,
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &A,
                            Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &M/*,
			    Dune::BlockVector<Dune::FieldVector<double,1> >& b*/
		    ){
                            //Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> &M_dtA_dune){
                            //Dune::DenseMatrix<Dune::FieldMatrix<double,1,1>> M_inverse) {

  using namespace Dune;


  auto gridView = basis->gridView();

  MatrixIndexSet occupationPattern;

  getOccupationPattern(basis, occupationPattern);

  occupationPattern.exportIdx(A);
  occupationPattern.exportIdx(M);
  //occupationPattern.exportIdx(M_dtA_dune);

  A = 0;
  M = 0;
  //M_dtA_dune=0;


  auto localView = basis->localView();
  auto localIndexSet = basis->localIndexSet();



  for (const auto &element : elements(gridView)) {

    localView.bind(element);
    localIndexSet.bind(localView);
    Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_A;
    assembleElementA(localView, elementMatrix_A);
    Dune::Matrix <FieldMatrix<double, 1, 1>> elementMatrix_M;
    assembleElementM(localView, elementMatrix_M);   


    for (size_t i = 0; i < elementMatrix_A.N(); i++) {


      auto row = localIndexSet.index(i);

      for (size_t j = 0; j < elementMatrix_A.M(); j++) {


        auto col = localIndexSet.index(j);


        A[row][col] += elementMatrix_A[i][j];
        //M_dtA_dune[row][col] += elementMatrix_A[i][j];
        M[row][col] += elementMatrix_M[i][j];


      }     
    }


  }

#if 0
  b.resize(basis->size());
  b= 0.0;
  double alpha = 0.1;
  auto flux = [&alpha](double T, const auto& x) { return std::abs(x[1]-2.0)<=0.25 ? alpha*(5.0-T) : 0.0; };
  fastBoundary(flux, *basis, gridView, b);
#endif
  /*
  storeMatrixMarket(A, "lapl-matrix.mm");
  storeMatrixMarket(M, "mass-matrix.mm");
  {
    std::fstream bOut("b.dat");
    bOut << b;
  }
  throw std::runtime_error("wrote matrices");
  */
}
