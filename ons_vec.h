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
struct mat2x2{
    union{
        struct{float a11; float a12; float a21; float a22;};
        struct {vec2 u; vec2 v; };
    };
};
struct mat3x3{
    union{
        struct{float a11; float a12;float a13; float a21; float a22;float a23;float a31; float a32; float a33;};
        struct {vec3 u; vec3 v;vec3 w;};
        struct {float a; float b; float c; float d; float e; float f; float g; float h; float i;};
        
    };
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

constexpr mat2x2 create_mat22(float a, float b, float c, float d)
{
    return {a,b,c,d};
}
constexpr mat2x2 operator+(const mat2x2& a, const mat2x2& b)
{
    return {a.a11+b.a11,a.a12+b.a12,a.a21+b.a21,a.a22+b.a22}; 
}
constexpr mat2x2 operator-(const mat2x2& a, const mat2x2& b)
{
    return {a.a11-b.a11,a.a12-b.a12,a.a21-b.a21,a.a22-b.a22}; 
}
constexpr mat2x2 operator*(const mat2x2& a, float x)
{
    return {a.a11*x,a.a12*x,a.a21*x,a.a22*x}; 
}
constexpr mat2x2 operator*(float x, const mat2x2& a)
{
    return {a.a11*x,a.a12*x,a.a21*x,a.a22*x}; 
}
constexpr vec2 operator*(const mat2x2& a, const vec2& v)
{
    return {a.a11*v.x+a.a12*v.y, a.a21*v.x+a.a22*v.y}; 
}
constexpr mat2x2 operator*(const mat2x2& a, const mat2x2& b)
{
    return {a.a11*b.a11+a.a12*b.a21, a.a11*b.a12+a.a12*b.a22,a.a21*b.a11+a.a22*b.a21, a.a21*b.a12+a.a22*b.a22}; 
}
constexpr mat2x2 transpose(const mat2x2& m)
{
    return {m.a11,m.a21,m.a12,m.a22};
}
constexpr float det(const mat2x2& m)
{
    return m.a11*m.a22-m.a12*m.a21; 
}
constexpr float det22(const mat2x2& m)
{
    return m.a11*m.a22-m.a12*m.a21; 
}
constexpr mat2x2 inv(const mat2x2& m)
{
    const float det_inv{(1/det(m))};
    return {det_inv*m.a22,-det_inv*m.a12,-det_inv*det_inv*m.a21,m.a11};
}
inline constexpr mat2x2 unity2x2{1,0,0,1};
inline constexpr mat2x2 zero2x2{0,0,0,0};
inline constexpr mat2x2 one2x2{1,1,1,1};


constexpr mat3x3 create_mat33(float a, float b, float c, float d,float e, float f, float g, float h, float i)
{
    return {a,b,c,d,e,f,g,h,i};
}
constexpr mat3x3 operator+(const mat3x3& a, const mat3x3& b)
{
    return {a.a11+b.a11,a.a12+b.a12,a.a13+b.a13,a.a21+b.a21,a.a22+b.a22,a.a23+b.a23,a.a31+b.a31,a.a32+b.a32,a.a33+b.a33}; 
}
constexpr mat3x3 operator-(const mat3x3& a, const mat3x3& b)
{
    return {a.a11-b.a11,a.a12-b.a12,a.a13-b.a13,a.a21-b.a21,a.a22-b.a22,a.a23-b.a23,a.a31-b.a31,a.a32-b.a32,a.a33-b.a33}; 
}
constexpr mat3x3 operator*(const mat3x3& a, float x)
{
    return {a.a11*x,a.a12*x,a.a13*x,a.a21*x,a.a22*x,a.a23*x,a.a31*x,a.a32*x,a.a33*x}; 
}
constexpr mat3x3 operator*(float x, const mat3x3& a)
{
    return {a.a11*x,a.a12*x,a.a13*x,a.a21*x,a.a22*x,a.a23*x,a.a31*x,a.a32*x,a.a33*x}; 
}
constexpr mat3x3 operator*(const mat3x3& a, const vec3& v)
{
    return {a.a11*v.x+a.a12*v.y+a.a13*v.z,a.a21*v.x+a.a22*v.y+a.a23*v.z,a.a31*v.x+a.a32*v.y+a.a33*v.z}; 
}
constexpr mat3x3 operator*(const mat3x3& a, const mat3x3& b)
{
    return {
    a.a11*b.a11+a.a12*b.a21+a.a13*b.a31, a.a11*b.a12+a.a12*b.a22+a.a13*b.a32, a.a11*b.a13+a.a12*b.a23+a.a13*b.a33,
    a.a21*b.a11+a.a22*b.a21+a.a23*b.a31, a.a21*b.a12+a.a22*b.a22+a.a23*b.a32, a.a21*b.a13+a.a22*b.a23+a.a23*b.a33,
    a.a31*b.a11+a.a32*b.a21+a.a33*b.a31, a.a31*b.a12+a.a32*b.a22+a.a33*b.a32, a.a31*b.a13+a.a32*b.a23+a.a33*b.a33,
    };
}
constexpr mat3x3 transpose(const mat3x3& m)
{
    return {
        m.a11,m.a21,m.a31,
        m.a12,m.a22,m.a32,
        m.a13,m.a23,m.a33
    };
}
constexpr float det(const mat3x3& m)
{
    return m.a11*det22({m.a22,m.a23,m.a32,m.a33})-m.a12*det22({m.a21,m.a23,m.a31,m.a33})+m.a13*det22({m.a21,m.a22,m.a31,m.a33});
}
constexpr float det33(const mat3x3& m)
{
    return m.a11*det22({m.a22,m.a23,m.a32,m.a33})-m.a12*det22({m.a21,m.a23,m.a31,m.a33})+m.a13*det22({m.a21,m.a22,m.a31,m.a33});
}
constexpr mat3x3 inv(const mat3x3& m)
{
    return create_mat33(
    m.e*m.i-m.f*m.h,m.c*m.h-m.b*m.i,m.b*m.f-m.c*m.e,
    m.f*m.g-m.d*m.i,m.a*m.i-m.c*m.g,m.c*m.d-m.a*m.f,
    m.d*m.h-m.e*m.g,m.b*m.g-m.a*m.h,m.a*m.e-m.b*m.d
    )*(1/det33(m)); 
}
inline constexpr mat3x3 unity3x3{1,0,0,0,1,0,0,0,1};
inline constexpr mat3x3 zero3x3{0,0,0,0,0,0,0,0,0};
inline constexpr mat3x3 one3x3{1,1,1,1,1,1,1,1,1};







} 


#endif