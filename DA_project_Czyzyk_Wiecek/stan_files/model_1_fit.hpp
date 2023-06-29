
// Code generated by stanc v2.29.0
#include <stan/model/model_header.hpp>
namespace model_1_fit_model_namespace {

using stan::model::model_base_crtp;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 42> locations_array__ = 
{" (found before start of program)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 11, column 4 to column 51)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 12, column 4 to column 49)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 13, column 4 to column 59)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 14, column 4 to column 59)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 15, column 4 to column 61)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 19, column 4 to column 32)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 20, column 4 to column 19)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 22, column 4 to column 28)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 23, column 4 to column 28)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 24, column 4 to column 55)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 25, column 4 to column 55)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 26, column 4 to column 57)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 45, column 4 to column 50)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 46, column 4 to column 50)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 47, column 4 to column 23)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 52, column 8 to column 128)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 53, column 8 to column 128)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 54, column 8 to column 145)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 51, column 20 to line 57, column 5)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 51, column 4 to line 57, column 5)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 30, column 4 to column 44)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 31, column 4 to column 43)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 33, column 4 to column 51)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 34, column 4 to column 52)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 35, column 4 to column 52)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 38, column 8 to column 125)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 39, column 8 to column 125)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 37, column 20 to line 40, column 5)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 37, column 4 to line 40, column 5)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 2, column 4 to column 10)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 3, column 26 to column 27)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 3, column 4 to column 29)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 4, column 21 to column 22)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 4, column 4 to column 24)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 5, column 29 to column 30)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 5, column 4 to column 32)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 6, column 29 to column 30)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 6, column 4 to column 32)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 45, column 47 to column 48)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 46, column 47 to column 48)",
 " (in '/home/projekt/stan_files/model_1_fit.stan', line 47, column 12 to column 13)"};




class model_1_fit_model final : public model_base_crtp<model_1_fit_model> {

 private:
  int N;
  std::vector<int> rating_difference;
  std::vector<int> time_control;
  std::vector<int> mistakes_white_model;
  std::vector<int> mistakes_black_model; 
  
 
 public:
  ~model_1_fit_model() { }
  
  inline std::string model_name() const final { return "model_1_fit_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.29.0", "stancflags = "};
  }
  
  
  model_1_fit_model(stan::io::var_context& context__,
                    unsigned int random_seed__ = 0,
                    std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "model_1_fit_model_namespace::model_1_fit_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 30;
      context__.validate_dims("data initialization","N","int",
           std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      
      
      current_statement__ = 30;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 31;
      stan::math::validate_non_negative_index("rating_difference", "N", N);
      current_statement__ = 32;
      context__.validate_dims("data initialization","rating_difference",
          "int", std::vector<size_t>{static_cast<size_t>(N)});
      rating_difference = 
        std::vector<int>(N, std::numeric_limits<int>::min());
      
      
      current_statement__ = 32;
      rating_difference = context__.vals_i("rating_difference");
      current_statement__ = 33;
      stan::math::validate_non_negative_index("time_control", "N", N);
      current_statement__ = 34;
      context__.validate_dims("data initialization","time_control","int",
           std::vector<size_t>{static_cast<size_t>(N)});
      time_control = std::vector<int>(N, std::numeric_limits<int>::min());
      
      
      current_statement__ = 34;
      time_control = context__.vals_i("time_control");
      current_statement__ = 35;
      stan::math::validate_non_negative_index("mistakes_white_model", "N", N);
      current_statement__ = 36;
      context__.validate_dims("data initialization","mistakes_white_model",
          "int", std::vector<size_t>{static_cast<size_t>(N)});
      mistakes_white_model = 
        std::vector<int>(N, std::numeric_limits<int>::min());
      
      
      current_statement__ = 36;
      mistakes_white_model = context__.vals_i("mistakes_white_model");
      current_statement__ = 37;
      stan::math::validate_non_negative_index("mistakes_black_model", "N", N);
      current_statement__ = 38;
      context__.validate_dims("data initialization","mistakes_black_model",
          "int", std::vector<size_t>{static_cast<size_t>(N)});
      mistakes_black_model = 
        std::vector<int>(N, std::numeric_limits<int>::min());
      
      
      current_statement__ = 38;
      mistakes_black_model = context__.vals_i("mistakes_black_model");
      current_statement__ = 39;
      stan::math::validate_non_negative_index("mistakes_white_pred", "N", N);
      current_statement__ = 40;
      stan::math::validate_non_negative_index("mistakes_black_pred", "N", N);
      current_statement__ = 41;
      stan::math::validate_non_negative_index("log_lik", "N", N);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1 + 1 + 1 + 1;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "model_1_fit_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ gamma_white = DUMMY_VAR__;
      current_statement__ = 1;
      gamma_white = in__.template read_constrain_lub<local_scalar_t__, 
                      jacobian__>(-0.0016, -0.0014, lp__);
      local_scalar_t__ gamma_black = DUMMY_VAR__;
      current_statement__ = 2;
      gamma_black = in__.template read_constrain_lub<local_scalar_t__, 
                      jacobian__>(0.0014, 0.0016, lp__);
      local_scalar_t__ time_controll_coeff_blitz = DUMMY_VAR__;
      current_statement__ = 3;
      time_controll_coeff_blitz = in__.template read_constrain_lub<
                                    local_scalar_t__, jacobian__>(1.13, 1.14,
                                    lp__);
      local_scalar_t__ time_controll_coeff_rapid = DUMMY_VAR__;
      current_statement__ = 4;
      time_controll_coeff_rapid = in__.template read_constrain_lub<
                                    local_scalar_t__, jacobian__>(0.76, 0.78,
                                    lp__);
      local_scalar_t__ time_controll_coeff_classic = DUMMY_VAR__;
      current_statement__ = 5;
      time_controll_coeff_classic = in__.template read_constrain_lub<
                                      local_scalar_t__, jacobian__>(0.20,
                                      0.22, lp__);
      std::vector<local_scalar_t__> time_controll_coeff =
         std::vector<local_scalar_t__>(3, DUMMY_VAR__);
      std::vector<local_scalar_t__> gammas =
         std::vector<local_scalar_t__>(2, DUMMY_VAR__);
      current_statement__ = 8;
      stan::model::assign(gammas, gamma_white,
        "assigning variable gammas", stan::model::index_uni(1));
      current_statement__ = 9;
      stan::model::assign(gammas, gamma_black,
        "assigning variable gammas", stan::model::index_uni(2));
      current_statement__ = 10;
      stan::model::assign(time_controll_coeff, time_controll_coeff_blitz,
        "assigning variable time_controll_coeff", stan::model::index_uni(1));
      current_statement__ = 11;
      stan::model::assign(time_controll_coeff, time_controll_coeff_rapid,
        "assigning variable time_controll_coeff", stan::model::index_uni(2));
      current_statement__ = 12;
      stan::model::assign(time_controll_coeff, time_controll_coeff_classic,
        "assigning variable time_controll_coeff", stan::model::index_uni(3));
      {
        current_statement__ = 21;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(
            stan::model::rvalue(gammas, "gammas", stan::model::index_uni(1)),
            -0.00153102, 0.0001));
        current_statement__ = 22;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(
            stan::model::rvalue(gammas, "gammas", stan::model::index_uni(2)),
            0.00152555, 0.0001));
        current_statement__ = 23;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(
            stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
              stan::model::index_uni(1)), 1.13878, 0.01));
        current_statement__ = 24;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(
            stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
              stan::model::index_uni(2)), 0.773998, 0.01));
        current_statement__ = 25;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(
            stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
              stan::model::index_uni(3)), 0.216708, 0.01));
        current_statement__ = 29;
        for (int ind = 1; ind <= N; ++ind) {
          current_statement__ = 26;
          lp_accum__.add(
            stan::math::poisson_log_lpmf<propto__>(
              stan::model::rvalue(mistakes_white_model,
                "mistakes_white_model", stan::model::index_uni(ind)),
              (stan::model::rvalue(time_controll_coeff,
                 "time_controll_coeff",
                 stan::model::index_uni(stan::model::rvalue(time_control,
                                          "time_control",
                                          stan::model::index_uni(ind)))) +
                (stan::model::rvalue(gammas, "gammas",
                   stan::model::index_uni(1)) *
                  stan::model::rvalue(rating_difference, "rating_difference",
                    stan::model::index_uni(ind))))));
          current_statement__ = 27;
          lp_accum__.add(
            stan::math::poisson_log_lpmf<propto__>(
              stan::model::rvalue(mistakes_black_model,
                "mistakes_black_model", stan::model::index_uni(ind)),
              (stan::model::rvalue(time_controll_coeff,
                 "time_controll_coeff",
                 stan::model::index_uni(stan::model::rvalue(time_control,
                                          "time_control",
                                          stan::model::index_uni(ind)))) +
                (stan::model::rvalue(gammas, "gammas",
                   stan::model::index_uni(2)) *
                  stan::model::rvalue(rating_difference, "rating_difference",
                    stan::model::index_uni(ind))))));
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "model_1_fit_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double gamma_white = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      gamma_white = in__.template read_constrain_lub<local_scalar_t__, 
                      jacobian__>(-0.0016, -0.0014, lp__);
      double gamma_black = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      gamma_black = in__.template read_constrain_lub<local_scalar_t__, 
                      jacobian__>(0.0014, 0.0016, lp__);
      double time_controll_coeff_blitz =
         std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 3;
      time_controll_coeff_blitz = in__.template read_constrain_lub<
                                    local_scalar_t__, jacobian__>(1.13, 1.14,
                                    lp__);
      double time_controll_coeff_rapid =
         std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 4;
      time_controll_coeff_rapid = in__.template read_constrain_lub<
                                    local_scalar_t__, jacobian__>(0.76, 0.78,
                                    lp__);
      double time_controll_coeff_classic =
         std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 5;
      time_controll_coeff_classic = in__.template read_constrain_lub<
                                      local_scalar_t__, jacobian__>(0.20,
                                      0.22, lp__);
      std::vector<double> time_controll_coeff =
         std::vector<double>(3, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> gammas =
         std::vector<double>(2, std::numeric_limits<double>::quiet_NaN());
      out__.write(gamma_white);
      out__.write(gamma_black);
      out__.write(time_controll_coeff_blitz);
      out__.write(time_controll_coeff_rapid);
      out__.write(time_controll_coeff_classic);
      if (stan::math::logical_negation((stan::math::primitive_value(
            emit_transformed_parameters__) || stan::math::primitive_value(
            emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 8;
      stan::model::assign(gammas, gamma_white,
        "assigning variable gammas", stan::model::index_uni(1));
      current_statement__ = 9;
      stan::model::assign(gammas, gamma_black,
        "assigning variable gammas", stan::model::index_uni(2));
      current_statement__ = 10;
      stan::model::assign(time_controll_coeff, time_controll_coeff_blitz,
        "assigning variable time_controll_coeff", stan::model::index_uni(1));
      current_statement__ = 11;
      stan::model::assign(time_controll_coeff, time_controll_coeff_rapid,
        "assigning variable time_controll_coeff", stan::model::index_uni(2));
      current_statement__ = 12;
      stan::model::assign(time_controll_coeff, time_controll_coeff_classic,
        "assigning variable time_controll_coeff", stan::model::index_uni(3));
      if (emit_transformed_parameters__) {
        out__.write(time_controll_coeff);
        out__.write(gammas);
      } 
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      std::vector<int> mistakes_white_pred =
         std::vector<int>(N, std::numeric_limits<int>::min());
      std::vector<int> mistakes_black_pred =
         std::vector<int>(N, std::numeric_limits<int>::min());
      Eigen::Matrix<double, -1, 1> log_lik =
         Eigen::Matrix<double, -1, 1>::Constant(N,
           std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 20;
      for (int ind = 1; ind <= N; ++ind) {
        current_statement__ = 16;
        stan::model::assign(mistakes_white_pred,
          stan::math::poisson_log_rng(
            (stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
               stan::model::index_uni(stan::model::rvalue(time_control,
                                        "time_control",
                                        stan::model::index_uni(ind)))) +
              (stan::model::rvalue(gammas, "gammas",
                 stan::model::index_uni(1)) *
                stan::model::rvalue(rating_difference, "rating_difference",
                  stan::model::index_uni(ind)))), base_rng__),
          "assigning variable mistakes_white_pred", stan::model::index_uni(ind));
        current_statement__ = 17;
        stan::model::assign(mistakes_black_pred,
          stan::math::poisson_log_rng(
            (stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
               stan::model::index_uni(stan::model::rvalue(time_control,
                                        "time_control",
                                        stan::model::index_uni(ind)))) +
              (stan::model::rvalue(gammas, "gammas",
                 stan::model::index_uni(2)) *
                stan::model::rvalue(rating_difference, "rating_difference",
                  stan::model::index_uni(ind)))), base_rng__),
          "assigning variable mistakes_black_pred", stan::model::index_uni(ind));
        current_statement__ = 18;
        stan::model::assign(log_lik,
          stan::math::poisson_log_lpmf<false>(
            stan::model::rvalue(mistakes_white_model, "mistakes_white_model",
              stan::model::index_uni(ind)),
            (stan::model::rvalue(time_controll_coeff, "time_controll_coeff",
               stan::model::index_uni(stan::model::rvalue(time_control,
                                        "time_control",
                                        stan::model::index_uni(ind)))) +
              (stan::model::rvalue(gammas, "gammas",
                 stan::model::index_uni(1)) *
                stan::model::rvalue(rating_difference, "rating_difference",
                  stan::model::index_uni(ind))))),
          "assigning variable log_lik", stan::model::index_uni(ind));
      }
      current_statement__ = 13;
      stan::math::check_greater_or_equal(function__, "mistakes_white_pred",
                                            mistakes_white_pred, 0);
      current_statement__ = 13;
      stan::math::check_less_or_equal(function__, "mistakes_white_pred",
                                         mistakes_white_pred, 50);
      current_statement__ = 14;
      stan::math::check_greater_or_equal(function__, "mistakes_black_pred",
                                            mistakes_black_pred, 0);
      current_statement__ = 14;
      stan::math::check_less_or_equal(function__, "mistakes_black_pred",
                                         mistakes_black_pred, 50);
      out__.write(mistakes_white_pred);
      out__.write(mistakes_black_pred);
      out__.write(log_lik);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_std_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ gamma_white = DUMMY_VAR__;
      gamma_white = in__.read<local_scalar_t__>();
      out__.write_free_lub(-0.0016, -0.0014, gamma_white);
      local_scalar_t__ gamma_black = DUMMY_VAR__;
      gamma_black = in__.read<local_scalar_t__>();
      out__.write_free_lub(0.0014, 0.0016, gamma_black);
      local_scalar_t__ time_controll_coeff_blitz = DUMMY_VAR__;
      time_controll_coeff_blitz = in__.read<local_scalar_t__>();
      out__.write_free_lub(1.13, 1.14, time_controll_coeff_blitz);
      local_scalar_t__ time_controll_coeff_rapid = DUMMY_VAR__;
      time_controll_coeff_rapid = in__.read<local_scalar_t__>();
      out__.write_free_lub(0.76, 0.78, time_controll_coeff_rapid);
      local_scalar_t__ time_controll_coeff_classic = DUMMY_VAR__;
      time_controll_coeff_classic = in__.read<local_scalar_t__>();
      out__.write_free_lub(0.20, 0.22, time_controll_coeff_classic);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"gamma_white", "gamma_black",
      "time_controll_coeff_blitz", "time_controll_coeff_rapid",
      "time_controll_coeff_classic", "time_controll_coeff", "gammas",
      "mistakes_white_pred", "mistakes_black_pred", "log_lik"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{}, std::vector<size_t>{
      }, std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(3)},
      std::vector<size_t>{static_cast<size_t>(2)},
      std::vector<size_t>{static_cast<size_t>(N)},
      std::vector<size_t>{static_cast<size_t>(N)},
      std::vector<size_t>{static_cast<size_t>(N)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "gamma_white");
    param_names__.emplace_back(std::string() + "gamma_black");
    param_names__.emplace_back(std::string() + "time_controll_coeff_blitz");
    param_names__.emplace_back(std::string() + "time_controll_coeff_rapid");
    param_names__.emplace_back(std::string() + "time_controll_coeff_classic");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 3; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "time_controll_coeff" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "gammas" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mistakes_white_pred" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mistakes_black_pred" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "log_lik" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "gamma_white");
    param_names__.emplace_back(std::string() + "gamma_black");
    param_names__.emplace_back(std::string() + "time_controll_coeff_blitz");
    param_names__.emplace_back(std::string() + "time_controll_coeff_rapid");
    param_names__.emplace_back(std::string() + "time_controll_coeff_classic");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 3; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "time_controll_coeff" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "gammas" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mistakes_white_pred" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "mistakes_black_pred" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "log_lik" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"gamma_white\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gamma_black\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_blitz\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_rapid\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_classic\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(3) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"gammas\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"mistakes_white_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"int\"}},\"block\":\"generated_quantities\"},{\"name\":\"mistakes_black_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"int\"}},\"block\":\"generated_quantities\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"generated_quantities\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"gamma_white\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gamma_black\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_blitz\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_rapid\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff_classic\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"time_controll_coeff\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(3) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"gammas\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"mistakes_white_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"int\"}},\"block\":\"generated_quantities\"},{\"name\":\"mistakes_black_pred\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(N) + ",\"element_type\":{\"name\":\"int\"}},\"block\":\"generated_quantities\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(N) + "},\"block\":\"generated_quantities\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  ((((1 + 1) + 1) + 1) + 1);
      const size_t num_transformed = (3 + 2);
      const size_t num_gen_quantities = 
  ((N + N) + N);
      std::vector<double> vars_vec(num_params__
       + (emit_transformed_parameters * num_transformed)
       + (emit_generated_quantities * num_gen_quantities));
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        vars_vec.data(), vars_vec.size());
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  ((((1 + 1) + 1) + 1) + 1);
      const size_t num_transformed = (3 + 2);
      const size_t num_gen_quantities = 
  ((N + N) + N);
      vars.resize(num_params__
        + (emit_transformed_parameters * num_transformed)
        + (emit_generated_quantities * num_gen_quantities));
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 5> names__{"gamma_white",
      "gamma_black", "time_controll_coeff_blitz",
      "time_controll_coeff_rapid", "time_controll_coeff_classic"};
      const std::array<Eigen::Index, 5> constrain_param_sizes__{1, 1, 
       1, 1, 1};
      const auto num_constrained_params__ = std::accumulate(
        constrain_param_sizes__.begin(), constrain_param_sizes__.end(), 0);
    
     std::vector<double> params_r_flat__(num_constrained_params__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < constrain_param_sizes__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(num_params_r__);
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits() 
    
};
}
using stan_model = model_1_fit_model_namespace::model_1_fit_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return model_1_fit_model_namespace::profiles__;
}

#endif


