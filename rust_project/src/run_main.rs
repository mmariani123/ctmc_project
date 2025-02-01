//#[derive(Debug)]

use std::env;

fn main() {
	
	//From the "Create_Parameter_Set.m" file:
	par.N        = 2 					//% Number of different starting genotypes
	par.a        = 24/(20/60) 			//% failure rate of male gametocytes (per day)
	par.b        = 24/(25/60) 			//% failure rate of female gametocytes (per day)
	par.r        = 0.08 				//% fertilization of male and female gametes (per day)
	par.mu_z     = 1 					//% death rate of zygotes (per day)
	par.sigma_z  = 24/19 				//% transformation rate of zygotes (per day)
	par.sigma_e  = 0.6 					//% transformation rate of ookinetes (per day)
	par.mu_e     = 1.4 					//% death rate of ookinetes (per day)
	par.mu_o     = 0 					//% death rate of oocysts (per day) **changed from paper**
	par.n        = 3e3 					//% number of sporozoites per oocyst
	par.p        = 0.2 					//% proportion of sporozoites that make it to salivary gland
	par.alpha    = 0.39 				//% fraction of male gametes that are viable
	par.beta     = 1 //%0.96    		//% fraction of female gametes that are viable
	par.eta      = 1 					//% number of female gametes per female gametocyte
	par.rho      = 8 					//% number of male gametes per male gametocytes
	par.k        = 1/7
	par.t0       = 10
	par.max_bias = 0.1 //%0.5 //%0.1
	par.NumSim   = 1 //%10000

	let mut params = HashMap::new();

    params.insert("N", 2);
    params.insert("a", 23/(20/60));
    params.insert("b", 24/(25/60));
    params.insert("r", 0.08);
	params.insert("mu_z", 1); 				
	params.insert("sigma_z", 24/19); 			
	params.insert("sigma_e", 0.6); 				
	params.insert("mu_e", 1.4) 				
	params.insert("mu_o", 0) 				
	params.insert("n", 3e3) 				
	params.insert("p", 0.2) 				
	params.insert("alpha", 0.39) 			
	params.insert("beta", 1) //%0.96    	
	params.insert("eta", 1) 				
	params.insert("rho", 8) 				
	params.insert("k", 1/7)
	params.insert("t0", 10)
	params.insert("max_bias", 0.1) //%0.5 //%0.1
	params.insert("NumSim", 1) //%10000
	
	N    	 = 2          //% Number of different starting genotypes
	a        = 24/(20/60) //% failure rate of male gametocytes (per day)
	b        = 24/(25/60) //% failure rate of female gametocytes (per day)
	r        = 0.08       //% fertilization of male and female gametes (per day)
	mu_z     = 1          //% death rate of zygotes (per day)
	sigma_z  = 24/19      //% transformation rate of zygotes (per day)
	sigma_e  = 0.6        //% transformation rate of ookinetes (per day)
	mu_e     = 1.4        //% death rate of ookinetes (per day)
	mu_o     = 0          //% death rate of oocysts (per day) **changed from paper**
	n        = 3e3        //% number of sporozoites per oocyst
	p        = 0.2        //% proportion of sporozoites that make it to salivary gland
	alpha    = 0.39       //% fraction of male gametes that are viable
	beta     = 1 //%0.96  //% fraction of female gametes that are viable
	eta      = 1          //% number of female gametes per female gametocyte
	rho      = 8          //% number of male gametes per male gametocytes
	k        = 1/7
	t0       = 10
	max_bias = 0.1 //%0.5 %0.1
	NumSim   = 1 //%10000	

	let mut percent_male = [[0u8; 4]; 6];
				percent[0][1] = 42;
				
	println!("{} days", percent[0][1]);
	
	//Make 2d array
	percent_male = [[f32; 2]; 2] 
	
	percent_male = [[.25],
					[.25]]
	
	println!("{} row1", percent_male[0][1]);
	println!("{} row2", percent_male[1][1]);
	
	if   N==1
		let mut bias = 0
	elif N==2
		let mut bias = [0,max_bias] //array
	elif N==3
		let mut bias = [0,0.1,0.5] //array 
	
	//transpose bias
	bias = np.transpose(bias)
	// Create a 2D array in row-major order: the rows of our 2D array are contiguous,
	// and the columns are strided
	let input_array = vec![ 1, 2, 3,
							4, 5, 6];
 
	// Treat our 6-element array as a 2D 3x2 array, and transpose it to a 2x3 array
	let mut output_array = vec![0; 6];
	transpose::transpose(&input_array, &mut output_array, 3, 2);

	// The rows have become the columns, and the columns have become the rows
	let expected_array =  vec![ 1, 4,
							2, 5,
							3, 6];
	assert_eq!(output_array, expected_array);

	// If we transpose our data again, we should get our original data back.
	let mut final_array = vec![0; 6];
	transpose::transpose(&output_array, &mut final_array, 2, 3);
	assert_eq!(final_array, input_array);
	
	//#.* means matrix product, if you don't write . 
	//#will Matlab product the numbers on the same position.

	//#a_vec = repmat(a,N,1).*(1 + bias);
	a_vec = np.matmul(tile(a,(N,1)), (1 + bias))
    
	//#b_vec = repmat(b,N,1).*(1 + bias);
	b_vec = np.matmul(tile(b,(N,1)), (1 + bias))

	//#alpha_vec = repmat(alpha,N,1).*(1 - bias);
	alpha_vec = np.matmul(tile(alpha,(N,1),(1-bias))

	//#beta_vec = repmat(beta,N,1).*(1 - bias);
	beta_vec = np.matmul(tile(beta,(N,1),[1 - bias])

	//#mu_z_vec = repmat(mu_z,N,1).*(1 + bias); 
	mu_z_vec = np.matmul(tile(z,(N,1)),(1 + bias))

	//#% zygote mortality rate for zygotes whose parents are of the same genotype

	//#mu_e_vec = repmat(mu_e,N,1).*(1 + bias); 
	mu_e_vec = np.matmul(tile(mu_e,(N,1)),(1+bias)

	//#% ookinete mortality rate for ookinetes " ...

	//#mu_o_vec = repmat(mu_o,N,1).*(1 + bias); 
	mu_o_vec = np.matmul(mu_o,(N,1),(1+bias))
	
	matrix multiplaction in rust
	use ndarray::arr2;

	fn main() {
    let a = arr2(&[[1, 2, 3],
                   [4, 5, 6]]);

    let b = arr2(&[[6, 3],
                   [5, 2],
                   [4, 1]]);

    println!("{}", a.dot(&b));
	}
	
	//multiple a scalar with a vector with a matrix
	use ndarray::{arr1, arr2, Array1};

	fn main() {
		let scalar = 4;

		let vector = arr1(&[1, 2, 3]);

		let matrix = arr2(&[[4, 5, 6],
                        [7, 8, 9]]);

		let new_vector: Array1<_> = scalar * vector;
		println!("{}", new_vector);

		let new_matrix = matrix.dot(&new_vector);
		println!("{}", new_matrix);
	}
	
	
	
	
	
	
	
	
	
	
	
	
    // Takes a reference and returns Option<&V>
    match contacts.get(&"Daniel") {
        Some(&number) => println!("Calling Daniel: {}", call(number)),
        _ => println!("Don't have Daniel's number."),
    }
	
    println!("Hello, world!");

	match variable_expression{
		case1=> {
		statement1
		}
		case2=> {
		statement2
		}
		_=> {
		Default statement
		}
	}

    const NUMBER: i32 = 20;
    const FLOATING_NUMBER: f32 = 20.0;
    println!("{}", NUMBER);

    println!("{}", FLOATING_NUMBER);

	enum Status {
		ACTIVE,
		INACTIVE,
	}

    let active = Status::ACTIVE;
    let inactive = Status::INACTIVE;
    println!("{:?}", active);
    println!("{:?}", inactive);

	for variable in expression-iterator{
		//code statements
	}

	fn main() {
		display_message();
	}

	fn display_message() {
		println!("Welcome to function tutorials in Rust")
	}

    println!("Please enter your Name");

    let args: Vec<String> = env::args().collect();

    let name = &args[1];
    let filename = &args[2];

    println!("Name is {}", name);

    let mut number = 0;
    loop {
        number += 1;
        if number == 3 {
            println! ("print continue statement, skipped printing the number");
            continue;
        }
        println!("{}", number);
        if number == 5 {
            println! ("Matched and exit from the loop");
            break; // Exit this loop
        }
    }

	let (first, second, third) = (false, 25, "John");
    println!("first: {} | second : {} | third = {}", first, second, third);

}

