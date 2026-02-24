#ifndef VEC3_H
#define VEC3_H

#include <cmath>

struct Vec3 {
    double x, y, z;

    Vec3(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}

    Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vec3& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    double magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }

    double magnitude_squared() const {
        return x * x + y * y + z * z;
    }

    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 operator/(double scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    Vec3& normalize() {
        double mag = magnitude();
        if (mag > 1e-12) {
            *this *= (1.0 / mag);
        }
        return *this;
    }

    static double apply_pbc(double x, double box_size) {
        if (box_size <= 0) return x;
        return x - box_size * round(x / box_size);
    }

    static Vec3 apply_pbc_vector(const Vec3& v, double box_size) {
        if (box_size <= 0) return v;
        return Vec3(apply_pbc(v.x, box_size), apply_pbc(v.y, box_size), apply_pbc(v.z, box_size));
    }
};

#endif // VEC3_H