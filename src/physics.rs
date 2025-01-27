extern crate nalgebra as na;
use integrators::verlet_integrate;
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

    // Vertlet
    pub fn update_bodies(&mut self, dt: f64) {

        let grav_bodies = self.gravitational_bodies.clone();
        
        let mut new_pos = Vec::<Vector3::<f64>>::new();
        let mut new_vel = Vec::<Vector3::<f64>>::new();
        let mut new_g_pos = Vec::<Vector3::<f64>>::new();
        let mut new_g_vel = Vec::<Vector3::<f64>>::new();
        
        for body in &self.gravitational_bodies {

            let old_acc = self.calculate_gravitational_acceleration(body.get_position());
            // Calculate the new position
            let new_position = body.get_position() + body.get_velocity() * dt + 0.5 * old_acc * dt * dt;

            // Calculate the new velocity
            let new_velocity = body.get_velocity() + 0.5 * (old_acc + self.calculate_gravitational_acceleration(new_position)) * dt;

            // Update the body's position and velocity
            new_g_pos.push(new_position);
            new_g_vel.push(new_velocity);
        }

        for body in &self.ballistic_bodies {

            let old_acc = self.calculate_gravitational_acceleration(body.get_position());
            // Calculate the new position
            let new_position = body.get_position() + body.get_velocity() * dt + 0.5 * old_acc * dt * dt;

            // Calculate the new velocity
            let new_velocity = body.get_velocity() + 0.5 * (old_acc + self.calculate_gravitational_acceleration(new_position)) * dt;

            // Update the body's position and velocity
            new_pos.push(new_position);
            new_vel.push(new_velocity);
        }

        for (i, body) in self.gravitational_bodies.iter_mut().enumerate() {
            body.set_position(new_g_pos[i]);
            body.set_velocity(new_g_vel[i]);
        }

        for (i, body) in self.ballistic_bodies.iter_mut().enumerate() {
            body.set_position(new_pos[i]);
            body.set_velocity(new_vel[i]);
        }

    }

}  // end of PhysicsWorld