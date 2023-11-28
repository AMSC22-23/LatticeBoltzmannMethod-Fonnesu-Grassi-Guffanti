#ifndef HH_BOUNDARY
#define HH_BOUNDARY

#include <vector>
#include <functional>

//We can use bounceBack when we have to deal with obstacle and 3D

class BounceBack
{

    virtual ~BounceBack() = default;
};

//TODO : choose parameters and implementation

class NEBB
{
   virtual std::vector<double> calc_top_wall_bounce() const = 0;
   virtual std::vector<double> calc_left_wall_bounce() const = 0;
   virtual std::vector<double> calc_right_wall_bounce() const = 0;
   virtual std::vector<double> calc_bottom_wall_bounce() const = 0;
   virtual std::vector<double> calc_top_left_corner_bounce() const = 0;
   virtual std::vector<double> calc_top_right_corner_bounce() const = 0;
   virtual std::vector<double> calc_bottom_left_corner_bounce() const = 0;
   virtual std::vector<double> calc_bottom_right_corner_bounce() const = 0;

   virtual ~NEBB() = default;
};

#endif