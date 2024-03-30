pub mod ekf;

extern crate nalgebra as na;

pub fn xyz_to_quaternion(euler:na::Vector3<f64>)->na::Quaternion<f64>
{
    let cos_a = (euler.x / 4.0).cos();
    let cos_b = (euler.y / 4.0).cos();
    let cos_r = (euler.z / 4.0).cos();
    let sin_a = (euler.x / 4.0).sin();
    let sin_b = (euler.y / 4.0).sin();
    let sin_r = (euler.z / 4.0).sin();

    na::Quaternion::<f64>::new(
        cos_a*cos_b*cos_r - sin_a*sin_b*sin_r, 
        sin_a*cos_b*cos_r + cos_a*sin_b*sin_r, 
        cos_a*sin_b*cos_r - sin_a*cos_b*sin_r, 
    cos_a*cos_b*sin_r + sin_a*sin_b*cos_r)
}

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