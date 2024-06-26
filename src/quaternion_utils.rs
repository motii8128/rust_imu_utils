extern crate nalgebra as na;

pub fn xyz_to_quaternion(euler:na::Vector3<f64>)->na::Quaternion<f64>
{
    let cos_a = (euler.x / 2.0).cos();
    let cos_b = (euler.y / 2.0).cos();
    let cos_r = (euler.z / 2.0).cos();
    let sin_a = (euler.x / 2.0).sin();
    let sin_b = (euler.y / 2.0).sin();
    let sin_r = (euler.z / 2.0).sin();

    na::Quaternion::<f64>::new(
        cos_a*cos_b*cos_r - sin_a*sin_b*sin_r, 
        sin_a*cos_b*cos_r + cos_a*sin_b*sin_r, 
        cos_a*sin_b*cos_r - sin_a*cos_b*sin_r, 
    cos_a*cos_b*sin_r + sin_a*sin_b*cos_r)
}

pub fn quaternion_to_xyz(w:f64, x:f64, y:f64, z:f64)->na::Vector3<f64>
{
    let y = (2.0*x*z + 2.0*y*w).asin();

    let x = if y.cos() == 0.0
        {
            let above = 2.0*y*z + 2.0*x*w;
            let below = 2.0*w*w + 2.0*y*y -1.0; 

            (above/below).atan()
        }
        else
        {
            let above = 2.0*y*z - 2.0*x*w;
            let below = 2.0*w*w + 2.0*z*z - 1.0;

            (-1.0*(above/below)).atan()
        };

    let z = if y.cos() == 0.0
        {
            0.0
        }
        else
        {
            let above = 2.0*x*y - 2.0*z*w;
            let below = 2.0*w*w + 2.0*x*x - 1.0;

            (-1.0*(above/below)).atan()
        };


    na::Vector3::<f64>::new(x, y, z)
}

pub fn zyx_to_quaternion(x:f64, y:f64, z:f64)->na::Quaternion<f64>
{
    let sin_x = (x/2.0).sin();
    let cos_x = (x/2.0).cos();
    let sin_y = (y/2.0).sin();
    let cos_y = (y/2.0).cos();
    let sin_z = (z/2.0).sin();
    let cos_z = (z/2.0).cos();

    na::Quaternion::new(
        sin_x*cos_y*cos_z - cos_x*sin_y*sin_z,
        sin_x*cos_y*sin_z + cos_x*sin_y*cos_z,
        -1.0*sin_x*sin_y*cos_z + cos_x*cos_y*sin_z,
        sin_x*sin_y*sin_z + cos_x*cos_y*cos_z
    )
}