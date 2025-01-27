extern crate nalgebra as na;
use na::{Vector3};

use super::*;

pub struct BallisticBody {
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

impl BallisticBody {
    pub fn new(mass: f64, position: Vector3<f64>, velocity: Vector3<f64>) -> Self {
        BallisticBody { mass, position, velocity }
    }
}

impl PhysicsBody for BallisticBody {
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