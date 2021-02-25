#include <iostream>
#include <cmath>
#include <string>

using std::string;
using std::sin;
using std::cos;
using std::sqrt;

bool equal_range (double lvalue, double rvalue){
    //lack of precision in a and w calculations force less precision here in order to correctly identify critically damped cases.
    const double range = .000000001; //1*10^-9
    return (lvalue <= rvalue + range && lvalue >= rvalue - range);
}

void print_usage(){
    std::cout << "Usage: Calculate values for solution of equations describing RLC circuits\
        \nSyntax: rlc -[p/s] -[n/f] \"R\" \"C\" \"L\" \"Src\"\
        \nFlags: -[p/s] (Parallel / Series RLC)\
        \n       -[n/f] (Natural / Forced response)\
        \n       -? (Print this msg)\
        \nR,C,L (Resistance (Ohms), Capacitance (Farads), Inductance(Henries))\
        \nSrc (voltage or current source value in a forced response RLC circuit)"<< std::endl;
}

int convert_args(char ** arg, double * R, double * C, double * L, double * src){
    if(arg[2] == nullptr || arg[3] == nullptr || arg[4] == nullptr){
        std::cout << "Invalid arguments use \"-?\" flag for help" << std::endl;
        return 0;
    }
    *R = std::stod(string(arg[3]));
    *C = std::stod(string(arg[4]));
    *L = std::stod(string(arg[5]));
    if(arg[6] != nullptr)*src = std::stod(string(arg[6]));
    else *src = 0;
    return 1;
}

//mode == 1; for series RLC
//mode == 0; for parallel RLC
void get_frequencies(double * a, double * w, double R, double C, double L, int mode){
    if(mode % 2 > 0) *a = R / (2.0 * L);
    else *a = 1.0 / (2.0*R*C);
    *w = 1.0 / sqrt(L * C);
}

void get_overdamped_exp(double * s1, double * s2, double a, double w){
    *s1 = (-1*a) + sqrt(a*a - w*w);
    //covers case where a and w are considered equal, they differ by < 10^-12
    //but w is slightly greater then a
    if(std::isnan(*s1)) *s1 = -1*a;
    *s2 = (-1*a) - sqrt(a*a - w*w);
}

void get_underdamped_period(double * T, double a, double w){
    *T = sqrt(w*w - a*a);
}

//mode == 3; for series forced response
//mode == 2; for parallel forced response
//mode == 1; for series RLC, capacitance isn't used
//mode == 0; for parallel RLC, inductance isn't used
//src argument ignored for mode == 0 and mode == 1
void get_overdamped_constant(double *c1, double *c2, double V_init, double I_init, double R, double L, double C, double s1, double s2, double src, int mode){
    double delta_init = 0; //the initial change in either voltage or current
    if(mode == 0){
        delta_init = (I_init * R + V_init) * -1.0 / L;
        *c2 = (delta_init - I_init*s1) / (s2 - s1);
        *c1 = I_init - *c2;
    }else if(mode == 1){
        delta_init = (I_init + V_init/R) * -1.0 / C;
        *c2 = (delta_init - V_init*s1) / (s2 - s1);
        *c1 = V_init - *c2;
    }else if(mode == 2){
        *c2 = ((V_init/L) - (s1 * I_init) + (s1 * src)) / (s2 - s1);
        *c1 = I_init - *c2 - src;
    }else{//mode == 3
        *c2 = ((I_init/C) - (s1 * V_init) + (s1 * src)) / (s2 - s1);
        *c1 = V_init - *c2 - src;
    }
}

//mode == 3; for series forced response
//mode == 2; for parallel forced response
//mode == 1; for series RLC, capacitance isn't used
//mode == 0; for parallel RLC, inductance isn't used
//src argument ignored for mode == 0 and mode == 1
void get_critcaldamped_constant(double *c1, double *c2, double V_init, double I_init, double R, double L, double C, double s, double src, int mode){
    double delta_init;
    if(mode == 0){
        delta_init = (I_init * R + V_init) * -1.0 / L;
        *c2 = I_init;
        *c1 = delta_init - (*c2 * s);
    }else if(mode == 1){
        delta_init = (I_init + V_init/R) * -1.0 / C;
        *c2 = V_init;
        *c1 = delta_init - (*c2 * s);
    }else if(mode == 2){
        *c2 = I_init - src;
        *c1 = (V_init / L) - s * (*c2);
    }else{//mode == 3
        *c2 = V_init - src;
        *c1 = (I_init / C) - s * (*c2);
    }
}

//mode == 3; for series forced response
//mode == 2; for parallel forced response
//mode == 1; for series RLC, capacitance isn't used
//mode == 0; for parallel RLC, inductance isn't used
//src argument ignored for mode == 0 and mode == 1
//a is neper frequency
//w is sqrt(a*a - w*w)
void get_underdamped_constant(double *c1, double *c2, double V_init, double I_init, double R, double L, double C, double a, double w, double src, int mode){
    double delta_init;
    if(mode == 0){
        delta_init = (I_init * R + V_init) * -1.0 / L;
        *c1 = I_init;
        *c2 = (delta_init + *c1 *a) / w;
    }else if(mode == 1){
        delta_init = (I_init + V_init/R) * -1.0 / C;
        *c1 = V_init;
        *c2 = (delta_init + *c1 *a) / w;
    }else if(mode == 2){
        *c1 = I_init - src;
        *c2 = ((V_init / L) + (a * *c1)) / w;
    }else{//mode == 3
        *c1 = V_init - src;
        *c2 = ((I_init / C ) + (a * *c1)) / w;
    }
}

int ask_for_constants(double *initials){
    std::cout << std::endl << "Calculate solution constants?    [Y\\N]" << std::endl;
    //delay
    std::cin.peek();
    while(std::cin.eof()){
        std::cin.peek();
    }
    char choice[100];
    std::cin.getline(choice, 2);
    if(choice[0] == 'Y' || choice[0] == 'y'){
        std::cout << std::endl << "Enter space seperated initial voltage followed by initial current" << std::endl;
        //delay
        std::cin.peek();
        while(std::cin.eof()){
            std::cin.peek();
        }
        std::cin.getline(choice, 100);
        //parse input
        string choice_str(choice);
        initials[0] = std::stod(choice_str.substr(0, choice_str.find(' ')));
        initials[1] = std::stod(choice_str.substr(choice_str.find(' ') + 1, choice_str.length()));
        return 1;
    }else return 0;
}

//finds and displays values related to the solutions of a natural Response RLC
//mode == 1; series
//mode == 0; paralell
int natural_calculations(double R, double C, double L, double src, int mode){
    double s1 = 0, s2 = 0, a = 0, w = 0, c1 = 0, c2 = 0;
    get_frequencies(&a, &w, R, C, L, mode);
    std::cout << "Found Neper frequency of: " << a << std::endl;
    std::cout << "Found resonant radian frequency of: " << w << std::endl << std::endl;
    if(equal_range(a, w)){
        std::cout << "Critically Damped:" << std::endl;
        get_overdamped_exp(&s1, &s2, a, w);
        std::cout << "Found exponent\n    s: " << s1 << std::endl;
    }else if(a > w){
        std::cout << "Overdamped:" << std::endl;
        get_overdamped_exp(&s1, &s2, a, w);
        std::cout << "Found exponents\n    s1: " << s1 << "\n    s2: " << s2 << std::endl;
    }else{//a < w
        std::cout << "Underdamped:" << std::endl;
        get_underdamped_period(&s1, a, w);
        std::cout << "Found damping frequency\n    w: " << s1 << std::endl;
    }
    double initials[2];
    if(ask_for_constants(initials) == 0)return 0;
    if(equal_range(a, w)){
        get_critcaldamped_constant(&c1, &c2, initials[0], initials[1], R, L, C, s1, src, mode);
    }else if(a > w){
        get_overdamped_constant(&c1, &c2, initials[0], initials[1], R, L, C, s1, s2, src, mode);
    }else{ //a < w
        get_underdamped_constant(&c1, &c2, initials[0], initials[1], R, L, C, a, s1, src, mode);
    }
    std::cout << "Found C1: " << c1 << "  C2: " << c2 << std::endl;
    return 1;
}


int main(int argc, char* argv[]){
    double R = 0, C = 0, L = 0, src = 0;
    std::cout.precision(20);
    //parse first flag
    std::cout << std::endl;
    //catch no argument usages
    if(argv[1] == nullptr){
        std::cout << "Flag required, use \"-?\" flag for help" << std::endl;
    }
    //display help msg
    else if(argv[1][1] == '?'){
        print_usage();
    }
    //parallel
    else if(argv[1][1] == 'p'){
        int mode = 0;
        std::cout << "Calculating for ";
        if(argv[2][1] == 'n') std::cout << "natural response ";
        else if(argv[2][1] == 'f'){
            std::cout << "forced response ";
            mode = 2;
        }
        else return 0;
        std::cout << "parallel RLC" << std::endl;
        if(! convert_args(argv, &R, &C, &L, &src)) return 0; //return if bad args
        natural_calculations(R, C, L, src, mode); //mode is 0 or 2
    }
    //series
    else if(argv[1][1] == 's'){
        int mode = 1;
        std::cout << "Calculating for ";
        if(argv[2][1] == 'n') std::cout << "natural response ";
        else if(argv[2][1] == 'f'){
            std::cout << "forced response ";
            mode = 3;
        }
        else return 0;
        std::cout << "series RLC" << std::endl;
        if(! convert_args(argv, &R, &C, &L, &src)) return 0; //return if bad args
        natural_calculations(R, C, L, src, mode); //mode is 1 or 3
    }
    //catch invalid inputs
    else{
        std::cout << "Flag required, use \"-?\" flag for help" << std::endl;
    }

    return 0;
}