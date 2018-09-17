# Numlib API Documentation

## Numlib Namespace Reference

### Typedefs

* `using Index = std::ptrdiff_t` 
* `template <typename T> using Band_mat = Band_matrix<T>`
* `template <typename T> using Vec = Matrix<T, 1>`
* `template <typename T> using Mat = Matrix<T, 2>`
* `template <typename T> using Cube = Matrix<T, 3>`
* `template <typename T> using Hypercube = Matrix<T, 4>`
* `template <typename T, Uplo_scheme Uplo> using Symm_mat = Packed_matrix<T, Uplo>`
* `template <typename T> using Upper_triang_mat = Packed_matrix<T, upper_triang>`
* `template <typename T> using Lower_triang_mat = Packed_matrix<T, lower_triang>`
* `template <typename T> using Sp_vec = Sparse_vector<T>`
* `template <typename T> using Sp_mat = Sparse_matrix<T>`

### Enumerations

* `enum Uplo_scheme {upper_triang, lower_triang}`
 
### Numlib::Constants

#### Mathematical Constants
* constexpr double **e** = 2.7182818284590452354 
* constexpr double **pi** = 3.14159265358979323846 
* constexpr double **sqrt2** = 1.41421356237309504880;

#### Metric Prefixes

* constexpr double **yotta** = 1.0e+24 
* constexpr double **zetta** = 1.0e+21 
* constexpr double **exa** = 1.0e+18 
* constexpr double **peta** = 1.0e+15 
* constexpr double **tera** = 1.0e+12 
* constexpr double **giga** = 1.0e+9 
* constexpr double **mega** = 1.0e+6 
* constexpr double **kilo** = 1.0e+3 
* constexpr double **hecto** = 1.0e+2 
* constexpr double **deca** = 10.0 
* constexpr double **one** = 1.0 
* constexpr double **deci** = 1.0e-1 
* constexpr double **centi** = 1.0e-2 
* constexpr double **milli** = 1.0e-3 
* constexpr double **micro** = 1.0e-6 
* constexpr double **nano** = 1.0e-9 
* constexpr double **pico** = 1.0e-12 
* constexpr double **femto** = 1.0e-15 
* constexpr double **atto** = 1.0e-18 
* constexpr double **zepto** = 1.0e-21 
* constexpr double **yocto** = 1.0e-24 

#### Physical Constants

* constexpr double **m_u** = 1.66053904000e-27
* constexpr double **N_A** = 6.02214085700e+23
* constexpr double **a_0** = 0.529177210670
* constexpr double **k** = 1.38064852000e-23
* constexpr double **G_0** = 7.74809173100e-5
* constexpr double **eps_0** = 8.85418781700e-12
* constexpr double **m_e** = 9.10938356000e-31
* constexpr double **eV** = 1.60217662080e-19
* constexpr double **ec** = 1.60217662080e-19
* constexpr double **F** = 9.64853328900e+4
* constexpr double **alpha** = 7.29735256640e-3
* constexpr double **R** = 8.31445980000
* constexpr double **E_h** = 4.35974465000e-18
* constexpr double **c_0** = 2.99792458e+8
* constexpr double **mu_0** = 1.25663706140e-6
* constexpr double **phi_0** = 2.06783383100e-15
* constexpr double **m_p_m_e** = 1.83615267389e+3
* constexpr double **G** = 6.67408000000e-11
* constexpr double **h** = 6.62607004000e-34
* constexpr double **h_bar** = 1.05457180000e-34
* constexpr double **m_p** = 1.67262189800e-27
* constexpr double **R_inf** = 1.09737315685e+7
* constexpr double **std_atm** = 1.01325000000e+5
* constexpr double **sigma** = 5.67036700000e-8

#### Conversion Factors:

* constexpr double **cal_to_J** = 4.184
* constexpr double **icm_to_kJ** = 1.19626564e-02
* constexpr double **icm_to_K** = 100.0 * h * c_0 / k
* constexpr double **J_to_icm** = 100.0 * h * c_0
* constexpr double **au_to_cm** = a_0 * 1.0e-8
* constexpr double **au_to_icm** = E_h / (h * c_0 * 100.0)
* constexpr double **au_to_s** = h_bar / E_h
* constexpr double **au_to_K** = E_h / k
* constexpr double **au_to_kg** = m_e
* constexpr double **au_to_kgm2** = m_u * a_0 * a_0 * 1.0e-20
* constexpr double **GHz_to_K** = giga * 4.79924470000e-11

### Classes

* class **Matrix< T, N >**
* class **Matrix_ref<T, N >**
* class **Matrix_base< N >**
* struct **Matrix_slice< N >**
* struct **slice**
* class **Band_matrix< T >**
* class **Packed_matrix< T >** 
* class **Sparse_vector< T >**
* class **Sparse_matrix< T >**


## Matrix Classes

### Numlib::Matrix<T, N>

N-dimensional dense matrix class using row-major storage order. The matrix 
class provides support for indexing, slicing and basic arithmetic operations.

**Template parameters:**
* **T** - The element type stored by the matrix
* **N** - The matrix order
 
#### Public Types

* using **size_type** = `typename Matrix_base<T, N>::size_type`
* using **value_type** = `typename Matrix_base<T, N>::value_type`
* using **iterator** = `typename std::vector<T>::iterator`
* using **const_iterator** = typename `std::vector<T>::const_iterator`

#### Public Member Functions

##### Constructors

* **Matrix**() = default
* template<typename... Exts> **Matrix**(Exts... exts)

##### Move Semantics

* **Matrix**(**Matrix**&&) = default
* **Matrix**& operator=(**Matrix**&&) = default

##### Copy Semantics

* **Matrix**(const **Matrix**&) = default
* **Matrix**& operator=(const **Matrix**&) = default
* template < typename U > **Matrix**(const **Matrix_ref**< U, N >&)

##### Initialize and Assign from List

* **Matrix**(Matrix_initializer< T, N >) 
* **Matrix**& operator=(Matrix_initializer< T, N >) 

##### Flat Element Access

* T* **data**() 
* const T* **data**() const 
