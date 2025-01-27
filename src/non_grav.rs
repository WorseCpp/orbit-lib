extern crate nalgebra as na;
use na::{Vector3};

const G: f64 = 6.67430e-11; // Gravitational constant in m^3 kg^-1 s^-2

pub struct GravitatingBody {
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

impl GravitatingBody {
    pub fn new(mass: f64, position: Vector3<f64>, velocity: Vector3<f64>) -> Self {
        GravitatingBody { mass, position, velocity }
    }

    pub fn get_acc(&self, r: Vector3<f64>) -> Vector3<f64> {
        let direction = self.position - r;
        let distance_squared = direction.magnitude_squared();
        let distance = distance_squared.sqrt();
        let force_magnitude = G * self.mass / distance_squared;
        direction.normalize() * force_magnitude
    }

    pub fn get_potential(&self, r: Vector3<f64>) -> f64 {
        let distance = (self.position - r).magnitude();
        -G * self.mass / distance
    }

    pub fn get_kinetic(&self) -> f64 {
        0.5 * self.mass * self.velocity.magnitude_squared()
    }

    pub fn get_grav_potential(&self, r: Vector3<f64>) -> f64 {
        let distance = (self.position - r).magnitude();
        -G * self.mass / distance
    }
}

impl Body for GravitatingBody {
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

