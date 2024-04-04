pub mod ekf;
pub mod quaternion_utils;

extern crate nalgebra as na;


pub fn convert_to_vector(x:f64, y:f64, z:f64)->na::Vector3<f64>
{
    na::Vector3::<f64>::new(x, y, z)
}

pub fn remove_gravity_from_posture(linear_accel:na::Vector3<f64>,euler:na::Vector3<f64>, gravity:f64)->na::Vector3<f64>
{
    let rm  = na::Rotation3::from_euler_angles(euler.x, euler.y, euler.z);

    let g = na::Vector3::new(0.0, 0.0, gravity);
    let removed = rm * g;

    

    // (removed.x + linear_accel.x, removed.y + linear_accel.y, removed.z - linear_accel.z)

    na::Vector3::<f64>::new(
        removed.x + linear_accel.x, 
        removed.y + linear_accel.y, 
        removed.z - linear_accel.z)
}