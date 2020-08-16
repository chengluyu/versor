use std::f32::consts::PI;

const RADIANS: f32 = PI / 180f32;

const DEGREES: f32 = 180f32 / PI;

/// Represents a spherical coordinate.
pub type SphericalCoordinate = (f32, f32);

/// Represents a 3-dimensional vector.
pub struct Vec3(f32, f32, f32);

/// Represents a quaternion.
#[derive(Copy, Clone)]
pub struct Quaternion(f32, f32, f32, f32);

fn hypot(a: f32, b: f32, c: f32, d: f32) -> f32 {
    (a * a + b * b + c * c + d * d).sqrt()
}

impl Quaternion {
    /// Computes the quaternion from a cartesian coordinate.
    pub fn from_cartesian(Vec3(x, y, z): Vec3) -> Quaternion {
        Quaternion(0f32, z, -y, x)
    }

    /// Computes the unit quaternion from a Euler rotation angles `(λ, φ, γ)`.
    pub fn from_angles(l: f32, p: f32, g: f32) -> Quaternion {
        let l = l * RADIANS / 2f32;
        let p = p * RADIANS / 2f32;
        let g = g * RADIANS / 2f32;
        let sl = l.sin();
        let cl = l.cos();
        let sp = p.sin();
        let cp = p.cos();
        let sg = g.sin();
        let cg = g.cos();
        Quaternion(
            cl * cp * cg + sl * sp * sg,
            sl * cp * cg - cl * sp * sg,
            cl * sp * cg + sl * cp * sg,
            cl * cp * sg - sl * sp * cg,
        )
    }

    /// Computes the Euler rotation angles `(λ, φ, γ)`.
    pub fn to_angles(&self) -> (f32, f32, f32) {
        let Quaternion(a, b, c, d) = self;
        (
            (2f32 * (a * b + c * d)).atan2(1f32 - 2f32 * (b * b + c * c)) * DEGREES,
            (2f32 * (a * c - d * b)).min(1f32).max(-1f32).asin() * DEGREES,
            (2f32 * (a * d + b * c)).atan2(1f32 - 2f32 * (c * c + d * d)) * DEGREES,
        )
    }

    /// Computes the length of `self`.
    pub fn length(&self) -> f32 {
        hypot(self.0, self.1, self.2, self.3)
    }

    /// Computes the dot product of `self` and `other`.
    pub fn dot(&self, rhs: &Quaternion) -> f32 {
        self.0 * rhs.0 + self.1 * rhs.1 + self.2 * rhs.2 + self.3 * rhs.3
    }

    /// Element-wisely multiplies `self` and `other`.
    pub fn multiply(&self, rhs: &Quaternion) -> Quaternion {
        let Quaternion(a1, b1, c1, d1) = self;
        let Quaternion(a2, b2, c2, d2) = rhs;
        Quaternion(
            a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
        )
    }
}

/// The type for interpolation functions.
pub type InterpolationFn = Box<dyn Fn(f32) -> Quaternion>;

/// Interpolate linearly between two quaternions.
pub fn linear_interpolate(a: Quaternion, b: Quaternion) -> InterpolationFn {
    let Quaternion(a1, b1, c1, d1) = a;
    let Quaternion(mut a2, mut b2, mut c2, mut d2) = b;
    a2 -= a1;
    b2 -= b1;
    c2 -= c1;
    d2 -= d1;
    Box::new(move |t: f32| -> Quaternion {
        let mut x = Quaternion(a1 + a2 * t, b1 + b2 * t, c1 + c2 * t, d1 + d2 * t);
        let l = x.length();
        x.0 /= l;
        x.1 /= l;
        x.2 /= l;
        x.3 /= l;
        x
    })
}

/// Interpolate between two quaternions.
pub fn interpolate(a: Quaternion, b: Quaternion) -> InterpolationFn {
    let Quaternion(a1, b1, c1, d1) = a;
    let Quaternion(mut a2, mut b2, mut c2, mut d2) = b;
    let mut dot = a.dot(&b);
    if dot < 0f32 {
        a2 = -a2;
        b2 = -b2;
        c2 = -c2;
        d2 = -d2;
        dot = -dot;
    }
    if dot > 0.9995 {
        return linear_interpolate(a, b);
    }
    let theta0 = dot.min(1f32).max(-1f32).acos();
    a2 -= a1 * dot;
    b2 -= b1 * dot;
    c2 -= c1 * dot;
    d2 -= d1 * dot;
    let l = hypot(a2, b2, c2, d2);
    a2 /= l;
    b2 /= l;
    c2 /= l;
    d2 /= l;
    Box::new(move |t: f32| -> Quaternion {
        let theta = theta0 * t;
        let s = theta.sin();
        let c = theta.cos();
        Quaternion(a1 * c + a2 * s, b1 * c + b2 * s, c1 * c + c2 * s, d1 * c + d2 * s)
    })
}

/// Returns Cartesian coordinates `[x, y, z]` given spherical coordinates `[λ, φ]`.
pub fn cartesian(e: SphericalCoordinate) -> (f32, f32, f32) {
    let l = e.0 * RADIANS;
    let p = e.1 * RADIANS;
    let cp = p.cos();
    (cp * l.cos(), cp * l.sin(), p.sin())
}

/// Returns the quaternion to rotate between two cartesian points on the sphere.
/// The value `alpha`, ranging from 0 to 1, is for tweening.
/// See https://github.com/Fil/versor/issues/8 for more.
pub fn delta(v0: Vec3, v1: Vec3, alpha: f32) -> Quaternion {
    fn cross(v0: &Vec3, v1: &Vec3) -> Vec3 {
        Vec3(
            v0.1 * v1.2 - v0.2 * v1.1,
            v0.2 * v1.0 - v0.0 * v1.2,
            v0.0 * v1.1 - v0.1 * v1.0
        )
    }

    fn dot(v0: &Vec3, v1: &Vec3) -> f32 {
        v0.0 * v1.0 + v0.1 * v1.1 + v0.2 * v1.2
    }
    
    let w = cross(&v0, &v1);
    let l = dot(&w, &w).sqrt();
    if l.is_nan() {
        Quaternion(1f32, 0f32, 0f32, 0f32)
    } else {
        let t = alpha * dot(&v0, &v1).min(1f32).max(-1f32).acos();
        let s = t.sin();
        Quaternion(t.cos(), w.2 / l * s, -w.1 / l * s, w.0 / l * s)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
