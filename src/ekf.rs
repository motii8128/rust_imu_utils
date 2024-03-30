extern crate nalgebra as na;

pub struct Axis6EKF
{
    estimation:na::Vector3<f64>,
    cov:na::Matrix3<f64>,
    gyro_noise:na::Matrix3<f64>,
    accel_noise:na::Matrix2<f64>,
    kalman_gain:na::Matrix3x2<f64>,
}

impl Axis6EKF {
    pub fn new(delta_time:f64)->Axis6EKF
    {
        let es_x = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let q = na::Matrix3::<f64>::zeros();
        let cov = na::Matrix3::<f64>::new(
            0.0174*delta_time*delta_time, 0.0, 0.0,
            0.0, 0.0174*delta_time*delta_time, 0.0,
            0.0, 0.0, 0.0174*delta_time*delta_time);

        let r = na::Matrix2::<f64>::new(
            0.0, 0.0, 
            0.0, 0.0);

        let k = na::Matrix3x2::<f64>::zeros();

        Axis6EKF{estimation:es_x,cov:cov,gyro_noise:q,accel_noise:r, kalman_gain:k}
    }
    pub fn run_ekf(&mut self,gyro_velocity:na::Vector3<f64>,accel:na::Vector3<f64>,delta_time:f64)->na::Vector3<f64>
    {
        let u = na::Vector3::<f64>::new(
            gyro_velocity.x*delta_time, 
            gyro_velocity.y*delta_time, 
            gyro_velocity.z*delta_time);

        self.gyro_noise = na::Matrix3::<f64>::new(
            0.0174*delta_time*delta_time, 0.0, 0.0,
            0.0, 0.0174*delta_time*delta_time, 0.0,
            0.0, 0.0, 0.0174*delta_time*delta_time,
        );

        self.accel_noise = na::Matrix2::<f64>::new(
            delta_time*delta_time, 0.0,
            0.0, delta_time
        );

        let jacob_f = self.calc_jacob(u);

        let _ = self.predict_x(u);

        let _ = self.predict_cov(jacob_f);

        let z = calc_observation_model(accel.x, accel.y, accel.z);

        let y_res = self.update_residual(z);

        let s = self.update_s();

        self.kalman_gain = self.update_kalman_gain(s);

        let _ = self.update_x(y_res);
        let _ = self.update_cov();

        self.estimation
    }
    fn h(&mut self)->na::Matrix2x3<f64>
    {
        na::Matrix2x3::<f64>::new(
            1.0, 0.0, 0.0, 
            0.0, 1.0, 0.0)
    }
    fn calc_jacob(&mut self, input_matrix:na::Vector3<f64>)->na::Matrix3<f64>
    {
        let cos_roll = self.estimation.x.cos();
        let sin_roll = self.estimation.x.sin();
        let cos_pitch = self.estimation.y.cos();
        let sin_pitch = self.estimation.y.sin();

        let m_11 = 1.0 + input_matrix.y*((cos_roll*sin_pitch)/cos_pitch) - input_matrix.z * ((sin_roll*sin_pitch)/cos_pitch);
        let m_12 = input_matrix.y*(sin_roll/(cos_pitch*cos_pitch))+input_matrix.z*((cos_roll/(cos_pitch*cos_pitch)));
        let m_21 = -1.0*input_matrix.y*sin_roll - input_matrix.z*cos_roll;
        let m_31 = input_matrix.y*(cos_roll/cos_pitch) - input_matrix.z*(sin_roll/cos_pitch);
        let m_32 = input_matrix.y*((sin_roll*sin_pitch)/(cos_pitch*cos_pitch))+input_matrix.z*((cos_roll*sin_pitch)/(cos_pitch*cos_pitch));

        na::Matrix3::<f64>::new(
            m_11, m_12, 0.0, 
            m_21, 1.0, 0.0, 
            m_31, m_32, 0.0)
    }
    fn predict_x(&mut self, input_matrix:na::Vector3<f64>)
    {
        let cos_roll = self.estimation.x.cos();
        let sin_roll = self.estimation.x.sin();
        let cos_pitch = self.estimation.y.cos();
        let sin_pitch = self.estimation.y.sin(); 


        self.estimation.x = self.estimation.x + input_matrix.x + input_matrix.y*((sin_roll*sin_pitch)/cos_pitch)+input_matrix.z*((cos_roll*sin_pitch)/cos_pitch);
        self.estimation.y = self.estimation.y + input_matrix.y * cos_roll - input_matrix.z*sin_roll;
        self.estimation.z = self.estimation.z + input_matrix.z + input_matrix.y*(sin_roll/cos_pitch) + input_matrix.z*(cos_roll/cos_pitch);
    }

    fn predict_cov(&mut self, jacob:na::Matrix3<f64>)
    {
        let t_jacob = jacob.transpose();
        self.cov = jacob*self.cov*t_jacob + self.gyro_noise;
    }
    fn update_residual(&mut self, obs:na::Vector2<f64>)->na::Vector2<f64>
    {
        let res = obs - self.h() * self.estimation;

        res
    }
    fn update_s(&mut self)->na::Matrix2<f64>
    {
        let h = self.h();
        let trans_h = h.transpose();

        h * self.cov * trans_h + self.accel_noise
    }
    fn update_kalman_gain(&mut self,s:na::Matrix2<f64>)->na::Matrix3x2<f64>
    {
        let h = self.h();
        let t_h = h.transpose();
        let inv_s = s.try_inverse().unwrap();
        let k = self.cov * t_h * inv_s;

        k
    }

    fn update_x(&mut self, residual:na::Vector2<f64>)
    {
        self.estimation = self.estimation + self.kalman_gain * residual;
    }
    fn update_cov(&mut self)
    {
        let i = na::Matrix3::<f64>::identity();
        let h = self.h();

        self.cov = (i - self.kalman_gain*h) * self.cov
    }
}



fn calc_observation_model(accel_x:f64, accel_y:f64, accel_z:f64)->na::Vector2<f64>
    {
        let x_ = match accel_z == 0.0 {
            true=>{
                if accel_y > 0.0
                {
                    std::f64::consts::PI / 2.0
                }
                else
                {
                    -1.0*std::f64::consts::PI / 2.0
                }
            },
            false=>(accel_y / accel_z).atan()
        };
        let y_ = match (accel_y*accel_y+accel_z*accel_z).sqrt() == 0.0 {
            true=>{
                if (-1.0*accel_x) > 0.0
                {
                    std::f64::consts::PI / 2.0
                }
                else
                {
                    -1.0*std::f64::consts::PI / 2.0
                }
            },
            false=>(-1.0*accel_x) / ((accel_y*accel_y+accel_z*accel_z).sqrt()).atan()
        };

        na::Vector2::<f64>::new(
            x_,
            y_
        )
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

pub struct Axis9EKF
{
    estimation:na::Vector3<f64>,
    cov:na::Matrix3<f64>,
    gyro_noise:na::Matrix3<f64>,
    accel_noise:na::Matrix3<f64>,
    kalman_gain:na::Matrix3<f64>,
}

impl Axis9EKF {
    pub fn new(delta_time:f64)->Axis9EKF
    {
        let es_x = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let q = na::Matrix3::<f64>::zeros();
        let cov = na::Matrix3::<f64>::new(
            0.0174*delta_time*delta_time, 0.0, 0.0,
            0.0, 0.0174*delta_time*delta_time, 0.0,
            0.0, 0.0, 0.0174*delta_time*delta_time);

        let r = na::Matrix3::<f64>::new(
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0);

        let k = na::Matrix3::<f64>::zeros();

        Axis9EKF{estimation:es_x,cov:cov,gyro_noise:q,accel_noise:r, kalman_gain:k}
    }
    pub fn run_ekf(&mut self,gyro_velocity:na::Vector3<f64>,accel:na::Vector3<f64>,magnet:na::Vector3<f64>,delta_time:f64)->na::Vector3<f64>
    {
        let u = na::Vector3::<f64>::new(
            gyro_velocity.x*delta_time, 
            gyro_velocity.y*delta_time, 
            gyro_velocity.z*delta_time);

        self.gyro_noise = na::Matrix3::<f64>::new(
            0.0174*delta_time*delta_time, 0.0, 0.0,
            0.0, 0.0174*delta_time*delta_time, 0.0,
            0.0, 0.0, 0.0174*delta_time*delta_time,
        );

        self.accel_noise = na::Matrix3::<f64>::new(
            delta_time*delta_time, 0.0, 0.0,
            0.0, delta_time*delta_time, 0.0,
            0.0, 0.0, delta_time*delta_time
        );

        let jacob_f = self.calc_jacob(u);

        let _ = self.predict_x(u);

        let _ = self.predict_cov(jacob_f);

        let z = observation_model(accel.x, accel.y, accel.z, magnet.x, magnet.y, magnet.z);

        let y_res = self.update_residual(z);

        let s = self.update_s();

        self.kalman_gain = self.update_kalman_gain(s);

        let _ = self.update_x(y_res);
        let _ = self.update_cov();

        self.estimation
    }
    fn calc_jacob(&mut self, input_matrix:na::Vector3<f64>)->na::Matrix3<f64>
    {
        let cos_roll = self.estimation.x.cos();
        let sin_roll = self.estimation.x.sin();
        let cos_pitch = self.estimation.y.cos();
        let sin_pitch = self.estimation.y.sin();

        let m_11 = 1.0 + input_matrix.y*((cos_roll*sin_pitch)/cos_pitch) - input_matrix.z * ((sin_roll*sin_pitch)/cos_pitch);
        let m_12 = input_matrix.y*(sin_roll/(cos_pitch*cos_pitch))+input_matrix.z*((cos_roll/(cos_pitch*cos_pitch)));
        let m_21 = -1.0*input_matrix.y*sin_roll - input_matrix.z*cos_roll;
        let m_31 = input_matrix.y*(cos_roll/cos_pitch) - input_matrix.z*(sin_roll/cos_pitch);
        let m_32 = input_matrix.y*((sin_roll*sin_pitch)/(cos_pitch*cos_pitch))+input_matrix.z*((cos_roll*sin_pitch)/(cos_pitch*cos_pitch));

        na::Matrix3::<f64>::new(
            m_11, m_12, 0.0, 
            m_21, 1.0, 0.0, 
            m_31, m_32, 0.0)
    }
    fn predict_x(&mut self, input_matrix:na::Vector3<f64>)
    {
        let cos_roll = self.estimation.x.cos();
        let sin_roll = self.estimation.x.sin();
        let cos_pitch = self.estimation.y.cos();
        let sin_pitch = self.estimation.y.sin(); 


        self.estimation.x = self.estimation.x + input_matrix.x + input_matrix.y*((sin_roll*sin_pitch)/cos_pitch)+input_matrix.z*((cos_roll*sin_pitch)/cos_pitch);
        self.estimation.y = self.estimation.y + input_matrix.y * cos_roll - input_matrix.z*sin_roll;
        self.estimation.z = self.estimation.z + input_matrix.z + input_matrix.y*(sin_roll/cos_pitch) + input_matrix.z*(cos_roll/cos_pitch);
    }

    fn predict_cov(&mut self, jacob:na::Matrix3<f64>)
    {
        let t_jacob = jacob.transpose();
        self.cov = jacob*self.cov*t_jacob + self.gyro_noise;
    }
    fn update_residual(&mut self, obs:na::Vector3<f64>)->na::Vector3<f64>
    {
        let res = obs - self.estimation;

        res
    }
    fn update_s(&mut self)->na::Matrix3<f64>
    {
        self.cov + self.accel_noise
    }
    fn update_kalman_gain(&mut self,s:na::Matrix3<f64>)->na::Matrix3<f64>
    {
        let inv_s = s.try_inverse().unwrap();
        let k = self.cov * inv_s;

        k
    }

    fn update_x(&mut self, residual:na::Vector3<f64>)
    {
        self.estimation = self.estimation + self.kalman_gain * residual;
    }
    fn update_cov(&mut self)
    {
        let i = na::Matrix3::<f64>::identity();

        self.cov = (i - self.kalman_gain) * self.cov
    }
}



fn observation_model(accel_x:f64, accel_y:f64, accel_z:f64, mag_x:f64, mag_y:f64, mag_z:f64)->na::Vector3<f64>
    {
        let x_ = match accel_z == 0.0 {
            true=>{
                if accel_y > 0.0
                {
                    std::f64::consts::PI / 2.0
                }
                else
                {
                    -1.0*std::f64::consts::PI / 2.0
                }
            },
            false=>(accel_y / accel_z).atan()
        };
        let y_ = match (accel_y*accel_y+accel_z*accel_z).sqrt() == 0.0 {
            true=>{
                if (-1.0*accel_x) > 0.0
                {
                    std::f64::consts::PI / 2.0
                }
                else
                {
                    -1.0*std::f64::consts::PI / 2.0
                }
            },
            false=>(-1.0*accel_x) / ((accel_y*accel_y+accel_z*accel_z).sqrt()).atan()
        };

        let cos_x = x_.cos();
        let sin_x = x_.sin();
        let cos_y = y_.cos();
        let sin_y= y_.sin();

        let z_ = match (mag_y*cos_x - mag_z*sin_x) == 0.0 {
            true=>{
                if (mag_x*cos_y + mag_y*sin_y*sin_x + mag_z*sin_y*cos_x) > 0.0
                {
                    std::f64::consts::PI / 2.0
                }
                else
                {
                    -1.0*std::f64::consts::PI / 2.0
                }
            },
            false=>((mag_x*cos_y + mag_y*sin_y*sin_x + mag_z*sin_y*cos_x) / (mag_y*cos_x - mag_z*sin_x)).atan()
        };

        na::Vector3::<f64>::new(
            x_,
            y_,
            z_
        )
}

// mag_x*cos_y + mag_y*sin_y*sin_x + mag_z*sin_y*cos_x