fn main() {
    let argv: Vec<String> = std::env::args().collect();
    std::process::exit(bwa_mem2_rs::generated::main_cpp::main(&argv));
}
