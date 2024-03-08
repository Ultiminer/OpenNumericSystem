#ifndef ONS_VEC_H_
#define ONS_VEC_H_
#include "quick_math.h"

namespace ONS
{
    
struct vec2{
    float x; float y; 
};
struct vec3{
    float x; float y; float z; 
};
struct vec2c
{
    QM::Complex x; QM::Complex y;
};
struct vec3c
{
    QM::Complex x; QM::Complex y;QM::Complex z;
};

constexpr vec2 operator+(const vec2& a, const vec2& b)
{
    return {a.x+b.x,a.y+b.y};
}
constexpr vec2 operator-(const vec2& a, const vec2& b)
{
    return {a.x-b.x,a.y-b.y};
}
constexpr vec2 operator*(float x, const vec2& a)
{
    return{x*a.x,x*a.y};
}
constexpr vec2 operator*( const vec2& a, float x)
{
    return{x*a.x,x*a.y};
}
constexpr float operator*(const vec2& a, const vec2& b)
{
    return a.x*b.x+a.y*b.y;
}
constexpr float cross(const vec2& a, const vec2& b)
{
    return a.x*b.y-a.y*b.x;
}
constexpr vec3 operator+(const vec3& a, const vec3& b)
{
    return {a.x+b.x,a.y+b.y,a.z+b.z};
}
constexpr vec3 operator-(const vec3& a, const vec3& b)
{
    return {a.x-b.x,a.y-b.y,a.z-b.z};
}
constexpr vec3 operator*(const vec3& a, float x)
{
    return {a.x*x,a.y*x,a.z*x};
}
constexpr vec3 operator*(float x, const vec3& a)
{
    return {a.x*x,a.y*x,a.z*x};
}
constexpr float operator*(const vec3& a, const vec3& b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
constexpr vec3 cross(const vec3& a, const vec3& b)
{
    return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}

constexpr vec2c operator+(const vec2c& a, const vec2c& b)
{
    return {a.x+b.x,a.y+b.y};
}
constexpr vec2c operator-(const vec2c& a, const vec2c& b)
{
    return {a.x-b.x,a.y-b.y};
}
constexpr vec2c operator*(const vec2c& a, float x)
{
    return {a.x*x,a.y*x};
}
constexpr vec2c operator*( float x,const vec2c& a)
{
    return {a.x*x,a.y*x};
}
constexpr QM::Complex operator*(const vec2c& a,const vec2c& b)
{
    return a.x.conjugate()*b.x+a.y.conjugate()*b.y; 
}
constexpr vec3c operator+(const vec3c& a, const vec3c& b)
{
    return {a.x+b.x,a.y+b.y,a.z+b.z};
}
constexpr vec3c operator-(const vec3c& a, const vec3c& b)
{
    return {a.x-b.x,a.y-b.y,a.z-b.z};
}
constexpr vec3c operator*(const vec3c& a, float x)
{
    return {a.x*x,a.y*x,a.z*x};
}
constexpr vec3c operator*(float x, const vec3c& a)
{
    return {a.x*x,a.y*x,a.z*x};
}
constexpr QM::Complex operator*(const vec3c& a,const vec3c& b)
{
    return a.x.conjugate()*b.x+a.y.conjugate()*b.y+a.z.conjugate()*b.z;
}






} 


#endif