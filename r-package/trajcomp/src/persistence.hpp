#include <Rcpp.h>
#include <cassert>
using namespace Rcpp;

namespace persistence{
  //basic vector math
  typedef std::vector<double> vec2;
  typedef std::vector<vec2> Trajectory;
  
  double length_sqr(const vec2 & v){
    return v[0] * v[0] + v[1] * v[1];
  }
  
  double length(const vec2 & v){
    return sqrt(length_sqr(v));
  }
  
  double dot(const vec2 & v, const vec2 & w){
    return v[0] * w[0] + v[1] * w[1];
  }
  
  //positive if v is left of w negative else
  double left_right(const vec2 & v, const vec2 & w){
    return v[0] * w[1] - v[1] * w[0];
  }
  
  //1.0 if v is left of w -1 else
  double lr_sgn(const vec2 & v, const vec2 & w){
    return (left_right(v, w) < 0) ? -1.0: 1.0;
  }
  
  double to_deg(double rad){
    return (rad * 180.0) / M_PI;
  }
  
  vec2 mult(const vec2 & v, double x){
    return {v[0] * x, v[1] * x};
  }
  
  void normalize(vec2 & v){
    auto len = length(v);
    if(len != 0){
      v[0] *= 1/len;
      v[1] *= 1/len;
    }
  }
  
  vec2 diff(const vec2 & v, const vec2 & w){
    return {v[0] - w[0], v[1] - w[1]};
  }
  
  vec2 sum(const vec2 & v, const vec2 & w){
    return {v[0] + w[0], v[1] + w[1]};
  }
  
  double distance(const vec2 & v, const vec2 & w){
    return length(diff(v,w));
  }
  
  double minimum_distance(const vec2 &v, const vec2 &w, const vec2 &p) {
    // Return minimum distance between line segment vw and point p
    const double l2 = length_sqr(diff(v, w));  // i.e. |w-v|^2 -  avoid a sqrt
    if (l2 == 0.0) return distance(p, v);   // v == w case
    // Consider the line extending the segment, parameterized as v + t (w - v).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-v) . (w-v)] / |w-v|^2
    // We clamp t from [0,1] to handle points outside the segment vw.
    const double t = std::max(0.0, std::min(1.0, dot(diff(p , v), diff(w , v)) / l2));
    const vec2 projection = sum(v ,mult( diff(w , v), t));  // Projection falls on the segment
    return distance(p, projection);
  }
  
  struct CurveElement{
    double val;
    int index;
    int comp_index;
    vec2 vertex;
  };
  
  typedef std::vector<CurveElement> Curve;
  
  struct Extrema{
    std::vector<int> min, max, used;
    std::set<int> max_set;
  };
  
  struct Component{
    int left, right;
    int min, max;
    int index;
    bool finished;
  };
  
  struct Result{
    std::vector<Component> comps;
    Curve pruned;
  };
  
  Curve traj_to_curve(const std::vector<vec2> &t){
    auto res = Curve(t.size());
    
    for(auto i = 0; i < t.size(); ++i){
      //edge cases
      if(i == 0 || i == t.size()-1){
        res[i] = {0, i, -1, t[i]};
      }else{
        auto v = t[i];
        auto last = t[i-1];
        auto v1 = diff(last, v);
        auto next = t[i+1];
        auto v2 = diff(v, next);
        
        auto len1 = length(v1);
        v1[0] *= 1/len1;
        v1[1] *= 1/len1;
        
        auto len2 = length(v2);
        v2[0] *= 1/len2;
        v2[1] *= 1/len2;
        
        auto value = dot(v1, v2);
        auto deg = to_deg( acos(value) * lr_sgn(v1, v2) );
        if(len1 == 0 || len2 == 0){
          deg = 0;
        }
        res[i] = {deg, i, -1, t[i]};
        
        //std::cout << i<<": "<<v1[0]<<", "<<v1[1] <<" | "<<v2[0]<<", "<<v2[1] <<" | "<<deg<<std::endl;
      }
    }
    return res;
  };
  
  void recalc_angles_inplace(Curve &curve){
    for(auto i = 0; i < curve.size(); ++i){
      //edge cases
      if(i == 0 || i == curve.size()-1){
        curve[i].val = 0;
      }else{
        auto v = curve[i].vertex;
        auto last = curve[i-1].vertex;
        auto v1 = diff(last, v);
        auto next = curve[i+1].vertex;
        auto v2 = diff(v, next);
        
        auto len1 = length(v1);
        v1[0] *= 1/len1;
        v1[1] *= 1/len1;
        
        auto len2 = length(v2);
        v2[0] *= 1/len2;
        v2[1] *= 1/len2;
        
        
        auto value = dot(v1, v2);
        auto deg = to_deg( acos(value)*lr_sgn(v1, v2) );
        if(len1 == 0 || len2 == 0){
          deg = 0;
        }
        curve[i].val = deg;
        curve[i].comp_index = -1;
      }
    }
  }
  
  Extrema find_extrema(const Curve &curve){
    const int fl_un = 0, fl_min = 1, fl_max = 2;
    Extrema res;
    
    auto d_min = curve[0], d_max = curve[0];
    auto search = fl_un;
    auto i_un = 0;
    
    for(auto i = 1; i < curve.size(); ++i){
      auto x = curve[i];
      if(search == fl_un){
        if(x.val > d_min.val){
          search = fl_max;
          d_max = x;
          res.min.push_back(i_un);
        }
        if(x.val < d_min.val){
          search = fl_min;
          d_min = x;
          res.max.push_back(i_un);
        }
      } else if(search == fl_min){
        if(x.val > d_min.val){
          search = fl_max;
          d_max = x;
          res.min.push_back(i-1);
        }else{
          d_min = x;
        }
      } else if(search == fl_max){
        if(x.val < d_max.val){
          search = fl_min;
          d_min = x;
          res.max.push_back(i-1);
        }else{
          d_max = x;
        }
      }
    }
    
    //if no min or max found
    if(res.min.size() == 0 || res.max.size() == 0){
      res.min.push_back(0);
      res.max.push_back(curve.size()-1);
    }else{
      //add end, its either min or max
      if(res.min.back() > res.max.back()){
        res.max.push_back(curve.size()-1);
      }else{
        res.min.push_back(curve.size()-1);
      }
    }
    
    res.used = std::vector<int>(curve.size());
    std::fill(res.used.begin(), res.used.end(), -1);
    
    res.max_set = std::set<int>(res.max.begin(), res.max.end());
    
    return res;
  }
  
  std::vector<Component> setup_components(const Extrema &extrema, Curve & curve){
    auto components = std::vector<Component>();
    int index = 0;
    for(auto m:extrema.min){
      curve[m].comp_index = index;
      components.push_back({m, m, m, -1, index, false});
      ++index;
    }
    return components;
  }
  
  int grow_component(Component & active_comp, const Curve &curve){
    auto l = active_comp.left - 1;
    auto r = active_comp.right + 1;
    
    //decide grow direction
    int x = 0;
    if((curve[l].val < curve[r].val || r == curve.size())  && l != -1  ){
      active_comp.left = l;
      x = l;
    }else{
      active_comp.right = r;
      x = r;
    }
    
    if(x == -1){
      x = 0;
    }
    if(x == curve.size()){
      x = curve.size()-1;
    }
    return x;
  }
  
  bool is_max(int index , const Extrema & extrema){
    return extrema.max_set.count(index) > 0 ;
  }
  
  bool is_collision(int index, const Extrema & extrema){
    return extrema.used[index] != -1;
  }

  Curve betaPruning( const std::vector<Component> &components, const Extrema &extrema,Curve &curve, double beta){
    for(auto comp: components){
      if(comp.max >= 0){
        //prune everything below beta
        if( (curve[comp.max].val - curve[comp.min].val) >= beta){
          curve[comp.min].comp_index = 1;
          curve[comp.max].comp_index = 1;
        }else{
          curve[comp.min].comp_index = -1;
          curve[comp.max].comp_index = -1;
        }
      }
    }
    
    curve[0].comp_index = 1;
    curve.back().comp_index = 1;
    
    Curve result;
    for(auto c: curve){
      if(c.comp_index == 1){
        result.push_back(c);
      }
    }
    return result;
  }
  
  void merge_collision(std::vector<Component> &components, const Extrema &extrema, Curve &curve,
                              Component &active_comp, int grow_index){
      Component & c2 = components[extrema.used[grow_index]];
      auto m1 = active_comp.min;
      auto m2 = c2.min;
      if (c2.min > active_comp.min){
        m1 = c2.min;
        m2 = active_comp.min;
      }
      
      c2.left = std::min(active_comp.left, c2.left);
      c2.right = std::max(active_comp.right, c2.right);
      c2.min = m2;
      c2.max = -1;
      c2.finished = false;
      
      active_comp.left = std::min(grow_index, m1);
      active_comp.right = std::max(grow_index, m1);
      active_comp.min = m1;
      active_comp.max = grow_index;
      active_comp.finished = true;
      
      curve[grow_index].comp_index = active_comp.index;
  }
  
  int next_component(const Component & active_comp, const std::vector<Component> & components){
    for(auto i = 0; i < components.size(); ++i){
      auto index = (i+active_comp.index+1) % components.size();
      
      if(components[index].finished != true){
        return index;
      }
    }
    return -1;
  }
  
  void mark_max_used(Component &active_comp, Extrema &extrema, int which){
    extrema.used[which] = active_comp.index;
  }
  
  std::vector<Component> build_components(Curve & curve, Extrema &extrema, int debug_iter = -1){
    auto components = setup_components(extrema, curve);
    
    //start with first component
    int active_index = 0;
    
    bool finished = false; 
    
    int it = 0;
    if(debug_iter < 0){
      debug_iter = curve.size()*(int)components.size();
    }
    int max_iter = std::min((int)curve.size()*(int)components.size(), debug_iter);
    
    
    while (! finished && ++it < max_iter){
      auto &active_comp = components[active_index];
      //decide grow direction
      auto grow_index = grow_component(active_comp, curve);
      
      //reaching a maximum?
      if(is_max(grow_index, extrema)){
        active_comp.finished = true;
        active_comp.max = grow_index;
        
        //maximum is already in use?
        if(is_collision(grow_index, extrema)){
          merge_collision(components, extrema, curve, active_comp, grow_index);
        }
        mark_max_used(active_comp, extrema, grow_index);
      }
      
      if(active_comp.left == 0 && active_comp.right == (int)curve.size()-1){
        finished = true;
        if(extrema.max_set.count(active_comp.left) > 0){
          active_comp.max = active_comp.left;
        }
        else if(extrema.max_set.count(active_comp.right) > 0){
          active_comp.max = active_comp.right;
        }
      }
      
      //returns index of next component or -1 if all finished
      if( (active_index = next_component(active_comp, components)) < 0){
        finished = true;
      }
    }
    return components;
  }
  
  Curve Persistence(Curve &curve, double beta){
    auto extrema = find_extrema(curve);
    if(curve.size() <3 || extrema.min.size() == 0 || extrema.max.size() == 0){
      return curve;
    }
    auto components = build_components(curve, extrema);
    return betaPruning(components, extrema, curve, beta);
  }
  
  Curve prune_curve_scale(const Curve &curve, double scale){
    if(curve.size() < 2){
      return curve;
    }
    Curve res;
    res.reserve(curve.size());
    res.push_back(curve[0]);
    
    for(auto i = 1; i < curve.size()-1;  ++i){
      auto it = curve[i], it2 = curve[i+1];
      //check distance
      if(distance(it.vertex, it2.vertex) < scale){
        //pruning
        if(it.val < it2.val){
          //avoid double insertion
          if(it.index != res.back().index)
            res.push_back(it);
        }else {
          res.push_back(it2);
        }
      }else{
        //no pruning
        res.push_back(it);
      }
    }
    res.push_back(curve.back());
    return res;
  }
  
  Curve prune_curve_dist_to_segment(const Curve &curve, double epsilon){
    if(curve.size() < 2){
      return curve;
    }
    Curve res;
    res.reserve(curve.size());
    res.push_back(curve[0]);
    vec2 last_pos = curve[0].vertex;
    for(auto i = 1; i < curve.size()-1;  ++i){
      auto it = curve[i], it2 = curve[i+1];
      //check distance
      if(minimum_distance(last_pos,it2.vertex, it.vertex) < epsilon){
        last_pos = it2.vertex;
        res.push_back(it2);
        ++i;
      }else{
        //no pruning
        last_pos = it.vertex;
        res.push_back(it);
      }
    }
    res.push_back(curve.back());
    return res;
  }
 
  Curve PersistenceMRS(Curve &curve,  double beta, int levels){
    auto _curve = curve;
    
    for(int i = 0; i < levels; ++i){
      if(i != 0){
        recalc_angles_inplace(_curve);
      }
      auto p_result = Persistence(_curve, beta);
      _curve = persistence::prune_curve_scale(p_result, pow(2,i));
    }
    return _curve;
  }
  
  Curve PersistenceSDS(Curve &curve,  double beta, double epsilon, int iterations){
    auto p_result = Persistence(curve, beta);
    Curve &_curve = p_result;
    
    for(int i = 0; i < iterations; ++i){
      _curve = prune_curve_dist_to_segment(_curve, epsilon);
    }
    return _curve;
  }

}//ENDOF NAMESPACE persistence

// [[Rcpp::export]]
NumericVector PersistenceCurve(NumericMatrix T) {
  std::vector<std::vector<double>> trajectory,query;
  
  for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
    trajectory.push_back({T(i,0),T(i,1)});
  
  auto curve = persistence::traj_to_curve(trajectory);
  auto res = std::vector<double>();
  for(auto &c: curve){
    res.push_back(c.val);
  }
  return wrap(res);
}

// [[Rcpp::export]]
NumericMatrix PersistenceBeta(NumericMatrix T, NumericVector Beta = 0) {
  persistence::Trajectory trajectory;
  try{
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    auto curve = persistence::traj_to_curve(trajectory);
    auto p_result = persistence::Persistence(curve, Beta(0));
    
    NumericMatrix resultMatrix = NumericMatrix(p_result.size(), 2) ;
    for(int i = 0; i < p_result.size(); i++){
      NumericVector temp  =  wrap(p_result[i].vertex);
      resultMatrix(i,_) = temp;
    }
    return resultMatrix;
  }catch(int e){
    std::cout<<"ERROR: "<< e <<std::endl;
  }
  return NumericMatrix();
}

// [[Rcpp::export]]
NumericMatrix PersistenceMRS(NumericMatrix T, NumericVector Beta = 0, NumericVector Levels = 5) {
  persistence::Trajectory trajectory;
  try{
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    auto curve = persistence::traj_to_curve(trajectory);
    auto p_result = persistence::PersistenceMRS(curve, Beta(0), Levels(0));
    
    NumericMatrix resultMatrix = NumericMatrix(p_result.size(), 2) ;
    for(int i = 0; i < p_result.size(); i++){
      NumericVector temp  =  wrap(p_result[i].vertex);
      resultMatrix(i,_) = temp;
    }
    return resultMatrix;
  }catch(int e){
    std::cout<<"ERROR: "<< e <<std::endl;
  }
  return NumericMatrix();
}

// [[Rcpp::export]]
NumericMatrix PersistenceSDS(NumericMatrix T, NumericVector Beta, NumericVector Epsilon, NumericVector Iterations = 3) {
  persistence::Trajectory trajectory;
  
    for (size_t i=0; i < T.nrow(); i++)		//@todo: remove copy by an adapter class @Martin
      trajectory.push_back({T(i,0),T(i,1)});
    
    auto curve = persistence::traj_to_curve(trajectory);
    auto p_result = persistence::PersistenceSDS(curve, Beta(0), Epsilon(0), Iterations(0));
    
    NumericMatrix resultMatrix = NumericMatrix(p_result.size(), 2) ;
    
    for(int i = 0; i < p_result.size(); i++){

      NumericVector temp  =  wrap(p_result[i].vertex);
      resultMatrix(i,_) = temp;
    }
    return resultMatrix;
}

// [[Rcpp::export]]
NumericVector Persistence_TestBars(NumericVector T) {
  persistence::Curve curve;
  
  for (size_t i=0; i < T.size(); i++)		//@todo: remove copy by an adapter class @Martin
    curve.push_back({T[i], (int)i, -1, {0,0}});
    
  auto extrema = persistence::find_extrema(curve);
  auto comps = persistence::build_components(curve, extrema);
  
  std::vector<int> bars;
  for(auto b : comps){
    if(b.max >= 0){
      bars.push_back( b.min);
      bars.push_back( b.max);
    }
  }
  return wrap(bars);
}

// [[Rcpp::export]]
NumericVector Persistence_TestComps(NumericVector T,  NumericVector it = -1) {
  persistence::Curve curve;
  
  for (size_t i=0; i < T.size(); i++)		//@todo: remove copy by an adapter class @Martin
    curve.push_back({T[i], (int)i, -1, {0,0}});
  
  auto extrema = persistence::find_extrema(curve);
  auto comps = persistence::build_components(curve, extrema, it(0));
  
  std::vector<int> bars;
  for(auto b : comps){
    std::cout <<"comp["<<b.index <<"]: "<<b.finished << ", "<<b.min<<", "<<b.max<<std::endl;
    bars.push_back( b.left);
    bars.push_back( b.min);
    bars.push_back( b.right);
  }
  return wrap(bars);
}