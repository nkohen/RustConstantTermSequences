mod dfao;
mod laurent_poly;
mod mod_int;
mod sequences;

use crate::dfao::DFAO;
use crate::laurent_poly::LaurentPoly;
use crate::sequences::{constant_term, constant_term_slow};

fn main() {
    let p = 5;
    let trinomial: LaurentPoly = LaurentPoly::from_string("x + 1 + x^-1", p);
    for i in 0..=20 {
        let tri_pow = trinomial.pow(&i);
        let mot_i = tri_pow.constant_term() - tri_pow.get_coefficient(&2);
        print!("{mot_i}, ");
    }

    let trinomial_pow = trinomial.pow(&9);
    let reduced = trinomial_pow.lambda_reduce();
    println!("\n{trinomial_pow} reduces to {reduced}");

    let one_minus_square: LaurentPoly = LaurentPoly::from_string("1-x^2", p);
    for i in 0..=20 {
        let fast_value = constant_term(&trinomial, &one_minus_square, &i);
        let slow_value = constant_term_slow(&trinomial, &one_minus_square, &i);
        if fast_value != slow_value {
            println!("{i} -- {fast_value}!= {slow_value}");
        }
    }
    println!("DONE!");

    let dfao = DFAO::poly_auto(&trinomial, &one_minus_square);
    for i in 0..=20 {
        print!("{}, ", dfao.compute_ct(i));
    }
    println!("\n{}", dfao.serialize());
    println!("{}", dfao.to_graphviz());

    dfao.save_png("./target/graph.png");
}
