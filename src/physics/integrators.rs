extern crate nalgebra as na;
use na::{Vector3};

use super::*;

pub fn verlet_integrate<F, B>(body: &mut B, dt: f64, get_acceleration: F)
where
    F: Fn(&Vector3<f64>) -> Vector3<f64>,
    B: PhysicsBody,
{
    let old_acc = get_acceleration(&body.get_position());
    // Calculate the new position
    let new_position = body.get_position() + body.get_velocity() * dt + 0.5 * old_acc * dt * dt;

    // Calculate the new velocity
    let new_velocity = body.get_velocity() + 0.5 * (old_acc + get_acceleration(&new_position)) * dt;

    // Update the body's position and velocity
    body.set_position(new_position);
    body.set_velocity(new_velocity);
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::grav::*;
    
    #[test]
    fn test_verlet_integration_energy_conservation() {
        let mass_earth = 5.972e24;
        let mass_satellite = 1000.0;
        let dt = 10.0f64;

        let mut satellite = GravitatingBody::new(
            mass_satellite,
            Vector3::new(7000.0e3, 0.0, 0.0),
            Vector3::new(0.0, 7.12e3, 0.0),
        );

        let earth = GravitatingBody::new(
            mass_earth,
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        );

        let initial_energy = satellite.get_kinetic() + satellite.get_mass() * earth.get_grav_potential_field(satellite.get_position());

        for _ in 0..1000 {
            verlet_integrate(&mut satellite, dt, |pos| earth.get_grav_acc(*pos));
            //println!("Energy: {}", satellite.get_kinetic() + satellite.get_mass() * earth.get_grav_potential_field(satellite.get_position()));
        }

        let final_energy = satellite.get_kinetic() + satellite.get_mass() * earth.get_grav_potential_field(satellite.get_position());

        let energy_difference = (final_energy - initial_energy).abs();
        let tolerance = initial_energy.abs() * 1.0e-5;

        assert!(energy_difference < tolerance, "Energy conservation test failed: difference = {}, tolerance = {}", energy_difference, tolerance);
    }
}
