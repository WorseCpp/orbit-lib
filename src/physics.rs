extern crate nalgebra as na;
use na::{Vector3};

mod grav;
mod non_grav;
mod integrators;

pub trait PhysicsBody {
    fn get_mass(&self) -> f64;
    fn get_position(&self) -> Vector3<f64>;
    fn get_velocity(&self) -> Vector3<f64>;
    fn set_position(&mut self, position: Vector3<f64>);
    fn set_velocity(&mut self, velocity: Vector3<f64>);
}