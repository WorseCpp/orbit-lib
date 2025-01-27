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
    fn boxed_clone(&self) -> Box<dyn PhysicsBody>;
}

impl Clone for Box<dyn PhysicsBody> {
    fn clone(&self) -> Box<dyn PhysicsBody> {
        self.boxed_clone()
    }
}

pub struct PhysicsWorld {
    gravitational_bodies: Vec<Box<dyn PhysicsBody>>,
    ballistic_bodies: Vec<Box<dyn PhysicsBody>>,
}

impl PhysicsWorld {
    pub fn calculate_gravitational_acceleration(&self, point: Vector3<f64>) -> Vector3<f64> {
        let mut acceleration = Vector3::zeros();
        for body in &self.gravitational_bodies {
            let direction = body.get_position() - point;
            let distance_squared = direction.norm_squared();
            if distance_squared > 0.0 {
                let force_magnitude = body.get_mass() / distance_squared;
                acceleration += direction.normalize() * force_magnitude;
            }
        }
        acceleration
    }

    pub fn update_bodies(&mut self, dt: f64) {

        let mut new_gravs = self.gravitational_bodies.clone();
        let mut new_balls = self.ballistic_bodies.clone();

        let grav_bodies = self.gravitational_bodies.clone();
        for body in &mut self.gravitational_bodies {
            
            let acc_fn = |pos: &Vector3<f64>| {
                let grav_bodies = &grav_bodies;
                self.calculate_gravitational_acceleration(*pos)
            };

            integrators::verlet_integrate(body.as_mut(), dt, acc_fn);
        }

        let ball_bodies = self.ballistic_bodies.clone();
        for body in &mut self.ballistic_bodies {
            
            let acc_fn = |pos: &Vector3<f64>| {
                let ball_bodies = &ball_bodies;
                self.calculate_gravitational_acceleration(*pos)
            };

            integrators::verlet_integrate(body.as_mut(), dt, acc_fn);
        }
    }

}  // end of PhysicsWorld