#include "pfasst/quadrature/uniform.hpp"

#include <stdexcept>
#include <vector>


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision>
    Uniform<precision>::Uniform(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      if (this->num_nodes < 2) {
        throw std::invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
      }
      this->compute_nodes();
      this->compute_weights();
    }

    template<typename precision>
    bool Uniform<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool Uniform<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }

    template<typename precision>
    void Uniform<precision>::compute_nodes()
    {
      this->nodes = std::vector<precision>(this->num_nodes, precision(0.0));
      for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = precision(j) / (this->num_nodes - 1);
      }
    }
  }  // ::pfasst::quadrature
  
  namespace quadrature
  {
    template<typename precision>
    Uniform_Right<precision>::Uniform_Right(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      if (this->num_nodes < 2) {
        throw std::invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
      }
      this->compute_nodes();
      this->compute_weights();
    }

    template<typename precision>
    bool Uniform_Right<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool Uniform_Right<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }
    
    template<typename precision>
    string Uniform_Right<precision>::print_summary() const
    {
      string ret = "Uniform (right) on (0,1.0] with " + std::to_string(this->get_num_nodes()) + " nodes";
      ret += "[";
      for(auto& d : this->get_nodes())
        ret += " " + std::to_string(d);
      ret += "]\n";
      return ret;
    }
    
    template<typename precision>
    void Uniform_Right<precision>::compute_nodes()
    {
      this->nodes = std::vector<precision>(this->num_nodes, precision(0.0));
      for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = precision(j+1) / (this->num_nodes);
      }
    }
    
    ///do the same extension as in gauss radau
    template<typename precision>
    void Uniform_Right<precision>::compute_weights()
    {
      IQuadrature<precision>::compute_weights();
      const size_t n = this->get_num_nodes();
      Matrix<precision> swap = this->q_mat;
      this->q_mat = Matrix<precision>::Zero(n+1, n+1);
      this->q_mat.block(1, 1, n, n) << swap;
      
      this->q_vec.insert(this->q_vec.begin(), precision(0));
      
      this->b_mat = Matrix<precision>(1, n + 1);
      for (size_t i = 0; i < n + 1; i++) {
        this->b_mat(0, i) = this->q_vec[i];
      }
    }
  }  // ::pfasst::quadrature
}  // ::pfasst
