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
    lifecycle_noise:na::Matrix3<f64>,
    obs_noise:na::Matrix3<f64>,
    kalman_gain:na::Matrix3<f64>,
}

impl Axis9EKF {
    pub fn new(delta_time:f64)->Axis9EKF
    {
        let es_x = na::Vector3::<f64>::new(0.0, 0.0, 0.0);
        let q = na::Matrix3::<f64>::new(
            0.01, 0.0, 0.0,
            0.0, 0.01, 0.0,
            0.0, 0.0, 0.01
        );
        let cov_ = na::Matrix3::<f64>::new(
            0.0174*delta_time*delta_time, 0.0, 0.0,
            0.0, 0.0174*delta_time*delta_time, 0.0,
            0.0, 0.0, 0.0174*delta_time*delta_time);

        let r = na::Matrix3::<f64>::new(
            0.1, 0.0, 0.0,
            0.0, 0.1, 0.0,
            0.0, 0.0, 0.1);

        let k = na::Matrix3::<f64>::zeros();

        Axis9EKF{estimation:es_x,cov:cov_,lifecycle_noise:q,obs_noise:r, kalman_gain:k}
    }
    pub fn run_ekf(&mut self,angular_velocity:na::Vector3<f64>,accel:na::Vector3<f64>,magnet:na::Vector3<f64>,delta_time:f64)->na::Vector3<f64>
    {
        let _ = self.predict_x(angular_velocity, delta_time);
        let obs = observation_model(accel, magnet);

        let est_jacob = self.estimation_jacob(angular_velocity, delta_time);
        let obs_jacob = self.observation_jacob();

        let pre_cov = self.predict_cov(est_jacob);
        let obs_cov = self.obs_cov(pre_cov, obs_jacob);

        let _ = self.calc_kalman_gain(pre_cov, obs_cov, obs_jacob);

        let _ = self.update_x(obs, obs_jacob);
        let _ = self.update_cov(obs_jacob, pre_cov);

        self.estimation
    }
    
    fn predict_x(&mut self, angular_velocity:na::Vector3<f64>, delta_time:f64)
    {
        let tan_y = self.estimation.y.tan();
        let sin_x = self.estimation.x.sin();
        let cos_x = self.estimation.x.cos();
        let cos_y = self.estimation.y.cos();

        self.estimation.x = self.estimation.x + (angular_velocity.x + angular_velocity.y*tan_y*sin_x + angular_velocity.z*tan_y*cos_x)* delta_time;
        self.estimation.y = self.estimation.y + (angular_velocity.y*cos_x - angular_velocity.z*sin_x) * delta_time;
        self.estimation.z = self.estimation.z + (angular_velocity.y*(sin_x/cos_y)+angular_velocity.z*(cos_x/cos_y))*delta_time;
    }

    fn estimation_jacob(&mut self, angular_velocity:na::Vector3<f64>, delta_time:f64)->na::Matrix3<f64>
    {
        let tan_y = self.estimation.y.tan();
        let sin_x = self.estimation.x.sin();
        let sin_y = self.estimation.y.sin();
        let cos_x = self.estimation.x.cos();
        let cos_y = self.estimation.y.cos();

        let m11 = 1.0 + (angular_velocity.x + angular_velocity.y*tan_y*cos_x - angular_velocity.z*tan_y*sin_x)*delta_time;
        let m12 = (angular_velocity.y*(sin_x/(cos_y*cos_y)) + angular_velocity.z*(cos_x/(cos_y*cos_y)))*delta_time;
        let m13 = 0.0;

        let m21 = (-1.0*angular_velocity.y*sin_x - angular_velocity.z*cos_x)*delta_time;
        let m22 = 1.0;
        let m23 = 0.0;

        let m31 = (angular_velocity.y*(cos_x/cos_y) - angular_velocity.z*(sin_x/cos_y))*delta_time;
        let m32 = (angular_velocity.y*sin_x*(sin_y/(cos_y*cos_y)) + angular_velocity.z*cos_x*(sin_y/(cos_y*cos_y)))*delta_time;
        let m33 = 1.0;

        na::Matrix3::<f64>::new(
            m11, m12, m13, 
            m21, m22, m23, 
            m31, m32, m33)
    }

    fn observation_jacob(&mut self)->na::Matrix3<f64>
    {
        na::Matrix3::<f64>::new(
            1.0, 0.0, 0.0, 
            0.0, 1.0, 0.0, 
            0.0, 0.0, 1.0)
    }

    fn predict_cov(&mut self, est_jacob:na::Matrix3<f64>)->na::Matrix3<f64>
    {
        let transpose_jacob = est_jacob.transpose();
        let p_cov = est_jacob*self.cov*transpose_jacob + self.lifecycle_noise;

        p_cov
    }

    fn obs_cov(&mut self, predict_cov:na::Matrix3<f64>, obs_jacob:na::Matrix3<f64>)->na::Matrix3<f64>
    {
        let transpose_jacob = obs_jacob.transpose();

        let o_cov = obs_jacob*predict_cov*transpose_jacob + self.obs_noise;

        o_cov
    }

    fn calc_kalman_gain(&mut self, predict_cov:na::Matrix3<f64>, obs_cov:na::Matrix3<f64>, obs_jacob:na::Matrix3<f64>)
    {
        let inversed_obs_cov = obs_cov.try_inverse().unwrap();
        let transpose_obs_jacob = obs_jacob.transpose();
        self.kalman_gain = predict_cov*transpose_obs_jacob*inversed_obs_cov;
    }

    fn update_x(&mut self, obs:na::Vector3<f64>, obs_jacob:na::Matrix3<f64>)
    {
        self.estimation = self.estimation + self.kalman_gain*(obs - obs_jacob*self.estimation);
    }

    fn update_cov(&mut self, obs_jacob:na::Matrix3<f64>, predict_cov:na::Matrix3<f64>)
    {
        let i = na::Matrix3::<f64>::new(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        );
        self.cov = (i - self.kalman_gain*obs_jacob)*predict_cov
    }
}



fn observation_model(accel_:na::Vector3<f64>, mag_:na::Vector3<f64>)->na::Vector3<f64>
{
    let mut result = na::Vector3::<f64>::new(0.0, 0.0, 0.0);

    result.x = ((-1.0*accel_.y) / (-1.0*accel_.z)).atan();

    result.y = (accel_.x / (accel_.y.powi(2) + accel_.z.powi(2)).sqrt()).atan();

    let above = mag_.x*result.y.cos() + mag_.y*result.y.sin()*result.x.sin() + mag_.z*result.y.sin()*result.x.cos();
    let below = mag_.y*result.x.cos() - mag_.z*result.x.sin();

    result.z = (above / below).atan();

    result    
}

// mag_x*cos_y + mag_y*sin_y*sin_x + mag_z*sin_y*cos_x