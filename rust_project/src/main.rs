use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {

    /// N
    #[arg(short, long, default_value_t = 1)]
    large_n: u8,
	
	#[arg(short, long, default_value_t = 1)]
    a: u8,
	
	#[arg(short, long, default_value_t = 1)]
    b: u8,
	
	#[arg(short, long, default_value_t = 1)]
    r: u8,
	
	#[arg(short, long, default_value_t = 1)]
    z_mu: u8,
	
	#[arg(short, long, default_value_t = 1)]
    isigma_z: u8,
	
	#[arg(short, long, default_value_t = 1)]
    jsigma_e: u8,
	
	#[arg(short, long, default_value_t = 1)]
    ee_mu: u8,
	
	#[arg(short, long, default_value_t = 1)]
    o_mu: u8,
	
	#[arg(short, long, default_value_t = 1)]
    n_small: u8,
	
	#[arg(short, long, default_value_t = 1)]
    p: u8,
	
	#[arg(short, long, default_value_t = 1)]
    u_alpha: u8,
	
	#[arg(short, long, default_value_t = 1)]
    x_beta: u8,
	
	#[arg(short, long, default_value_t = 1)]
    y_eta: u8,
	
	#[arg(short, long, default_value_t = 1)]
    w_rho: u8,
	
	#[arg(short, long, default_value_t = 1)]
    current_time: u8,
	
	#[arg(short, long, default_value_t = 1)]
    time_zero: u8,
	
	#[arg(short, long, default_value_t = 1)]
    max_bias: u8,
	
	#[arg(short, long, default_value_t = 1)]
    sim_num: u8,
	
	//par.N        = 2 				
	//par.a        = 24/(20/60) 		
	//par.b        = 24/(25/60) 		
	//par.r        = 0.08 			
	//par.mu_z     = 1 				
	//par.sigma_z  = 24/19 			
	//par.sigma_e  = 0.6 				
	//par.mu_e     = 1.4 				
	//par.mu_o     = 0 				
	//par.n        = 3e3 				
	//par.p        = 0.2 				
	//par.alpha    = 0.39 			
	//par.beta     = 1 //%0.96    	
	//par.eta      = 1 				
	//par.rho      = 8 				
	//par.k        = 1/7
	//par.t0       = 10
	//par.max_bias = 0.1 //%0.5 //%0.1
	//par.NumSim   = 1 //%10000
	
	//let matches = App::new("clap")
    //.arg(Arg::with_name("Fmin")
    //   .required(false)
    //    .takes_value(true)
    //    .short("Fmin")
    //    .multiple(false)
    //    .possible_values(&["min"])
    //)
    //.get_matches();
	
}

fn main() {
    
	let args = Args::parse();

    for _ in 0..args.sim_num {
        println!("Number of simulations {}!", args.sim_num);
	}
	
}
