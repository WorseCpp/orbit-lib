extern crate nalgebra as na;
use na::{Vector3};

use super::*;

const G: f64 = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2

pub struct GravitatingBody {
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>
}

impl GravitatingBody {
    pub fn new(mass: f64, position: Vector3<f64>, velocity: Vector3<f64>) -> Self {
        GravitatingBody { mass, position, velocity}
    }

    pub fn get_grav_acc(&self, r: Vector3<f64>) -> Vector3<f64> {
        let direction = self.position - r;
        let distance_squared = direction.magnitude_squared();
        let distance = distance_squared.sqrt();
        let force_magnitude = G * self.mass / distance_squared;
        direction.normalize() * force_magnitude
    }

    pub fn get_grav_potential_field(&self, r: Vector3<f64>) -> f64 {
        let distance = (self.position - r).magnitude();
        -G * self.mass / distance
    }
    pub fn get_kinetic(&self) -> f64 {
        0.5 * self.mass * self.velocity.magnitude_squared()
    }
}

impl PhysicsBody for GravitatingBody {

    fn boxed_clone(&self) -> Box<dyn PhysicsBody> {
        Box::new(GravitatingBody {
            mass: self.mass,
            position: self.position,
            velocity: self.velocity,
        })
    }
    
    fn get_mass(&self) -> f64 {
        self.mass
    }

    fn get_position(&self) -> Vector3<f64> {
        self.position
    }

    fn get_velocity(&self) -> Vector3<f64> {
        self.velocity
    }

    fn set_position(&mut self, position: Vector3<f64>) {
        self.position = position;
    }

    fn set_velocity(&mut self, velocity: Vector3<f64>) {
        self.velocity = velocity;
    }


}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acceleration_at_surface_of_earth() {
        let earth_mass = 5.972e24; // in kg
        let earth_radius = 6.371e6; // in meters
        let earth_position = Vector3::new(0.0, 0.0, 0.0);
        let earth_velocity = Vector3::new(0.0, 0.0, 0.0);

        let earth = GravitatingBody::new(earth_mass, earth_position, earth_velocity);
        let surface_position = Vector3::new(earth_radius, 0.0, 0.0);

        let acceleration = earth.get_grav_acc(surface_position);
        let expected_acceleration = 9.81; // Approximate acceleration due to gravity at Earth's surface in m/s^2

        assert!((acceleration.magnitude() - expected_acceleration).abs() < 0.1);
    }
}

