extern crate nalgebra as na;

use na::{Vector3};
use grav::GravitatingBody;
use non_grav::BallisticBody;

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
    fn get_kinetic(&self) -> f64 {
        0.5 * self.get_mass() * self.get_velocity().norm_squared()
    }
}

impl Clone for Box<dyn PhysicsBody> {
    fn clone(&self) -> Box<dyn PhysicsBody> {
        self.boxed_clone()
    }
}

pub struct PhysicsWorld {
    gravitational_bodies: Vec<GravitatingBody>,
    ballistic_bodies: Vec<Box<dyn PhysicsBody>>,
}

impl PhysicsWorld {
    pub fn calculate_gravitational_acceleration(&self, point: Vector3<f64>) -> Vector3<f64> {
        let mut acceleration = Vector3::zeros();
        for body in &self.gravitational_bodies {
            acceleration += body.get_grav_acc(point);
        }
        acceleration
    }

    pub fn calculate_gravitational_acceleration_exclude(&self, point: Vector3<f64>, exclude: &GravitatingBody) -> Vector3<f64> {
        let mut acceleration = Vector3::zeros();
        for body in &self.gravitational_bodies {
            
            if std::ptr::eq(body, exclude) {
                continue;
            }
            
            acceleration += body.get_grav_acc(point);
        }
        acceleration
    }

    pub fn calculate_gravitational_potential_field(&self, point: Vector3<f64>) -> f64 {
        let mut potential_field = 0.0;
        for body in &self.gravitational_bodies {
            let distance = (body.get_position() - point).norm();
            if distance > 0.01 {
                potential_field += body.get_grav_potential_field(point);
            }
        }
        potential_field
    }
    pub fn calculate_gravitational_potential_field_exclude(&self, point: Vector3<f64>, exclude: &GravitatingBody) -> f64 {
        let mut potential_field = 0.0;
        for body in &self.gravitational_bodies {
            let distance = (body.get_position() - point).norm();
            if std::ptr::eq(body, exclude) {
                continue;
            }
            if distance > 0.01 {
                potential_field += body.get_grav_potential_field(point);
            }
        }
        potential_field
    }

    

    pub fn total_energy(&self) -> f64 {
        let mut total_energy = 0.0;

        for body in &self.gravitational_bodies {
            total_energy += body.get_kinetic();
            for other_body in &self.gravitational_bodies {
                if std::ptr::eq(body, other_body) {
                    break;
                }
                total_energy += body.get_grav_potential_field(other_body.get_position()) * other_body.get_mass();
            }
        }

        for body in &self.ballistic_bodies {
            total_energy += body.get_kinetic();
            total_energy += self.calculate_gravitational_potential_field(body.get_position()) * body.get_mass();
        }

        total_energy
    }

    // Vertlet
    pub fn update_bodies(&mut self, dt: f64) {

        
        let mut new_pos = Vec::<Vector3::<f64>>::new();
        let mut new_vel = Vec::<Vector3::<f64>>::new();
        let mut new_g_pos = Vec::<Vector3::<f64>>::new();
        let mut new_g_vel = Vec::<Vector3::<f64>>::new();
        
        for body in &self.gravitational_bodies {

            let old_acc = self.calculate_gravitational_acceleration_exclude(body.get_position(), body);
            // Calculate the new position
            let new_position = body.get_position() + body.get_velocity() * dt + 0.5 * old_acc * dt * dt;

            // Calculate the new velocity
            let new_velocity = body.get_velocity() + 0.5 * (old_acc + self.calculate_gravitational_acceleration_exclude(new_position, body)) * dt;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_physics_world_with_planets() {
        let sun_mass = 1.989e30;
        let sun_position = Vector3::new(0.0, 0.0, 0.0);
        let sun_velocity = Vector3::new(0.0, 0.0, 0.0);
        let sun = GravitatingBody::new(sun_mass, sun_position, sun_velocity);

        let earth_mass = 5.972e24;
        let earth_position = Vector3::new(1.496e11, 0.0, 0.0); // 1 AU
        let earth_velocity = Vector3::new(0.0, 29.78e3, 0.0); // Earth's orbital speed
        let earth = GravitatingBody::new(earth_mass, earth_position, earth_velocity);

        let mars_mass = 6.4171e23;
        let mars_position = Vector3::new(2.279e11, 0.0, 0.0); // 1.524 AU
        let mars_velocity = Vector3::new(0.0, 24.077e3, 0.0); // Mars' orbital speed
        let mars = GravitatingBody::new(mars_mass, mars_position, mars_velocity);

        let mut world = PhysicsWorld {
            gravitational_bodies: vec![sun, earth, mars],
            ballistic_bodies: vec![],
        };

        let initial_energy = world.total_energy();
        println!("Initial Energy: {}", initial_energy);
        let dt = 60.0 * 60.0 * 24.0; // 1 day in seconds
        let days_in_year = 365;
        for _ in 0..days_in_year {
            world.update_bodies(dt);

        }
        
        let final_energy = world.total_energy();
        assert!(((initial_energy / final_energy) -1.0 ).abs() < 0.001, "Energy is not conserved");
        
    }

}
