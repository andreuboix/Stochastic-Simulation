#include <iostream>
#include <queue>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
using P = pair<double,int>;
using PQ = priority_queue<P, vector<P>, greater<P>>;
// We use a priority_queue such that we use a min-heap (the .front() gives us the smallest value)
// This translates for us that we get the earliest event (the event with smallest time value)

// Approximate inverse CDF for the standard normal distribution (using Abramowitz and Stegun)
double normal_quantile(double p) {
    if (p <= 0.0 || p >= 1.0) {
        throw std::invalid_argument("Percentile must be between 0 and 1");
    }

    // Coefficients for the approximation
    static const double a1 = -3.969683028665376e+01;
    static const double a2 =  2.209460984245205e+02;
    static const double a3 = -2.759285104469687e+02;
    static const double a4 =  1.383577518672690e+02;
    static const double a5 = -3.066479806614716e+01;
    static const double a6 =  2.506628277459239e+00;

    static const double b1 = -5.447609879822406e+01;
    static const double b2 =  1.615858368580409e+02;
    static const double b3 = -1.556989798598866e+02;
    static const double b4 =  6.680131188771972e+01;
    static const double b5 = -1.328068155288572e+01;

    static const double c1 = -7.784894002430293e-03;
    static const double c2 = -3.223964580411365e-01;
    static const double c3 = -2.400758277161838e+00;
    static const double c4 = -2.549732539343734e+00;
    static const double c5 =  4.374664141464968e+00;
    static const double c6 =  2.938163982698783e+00;

    static const double d1 =  7.784695709041462e-03;
    static const double d2 =  3.224671290700398e-01;
    static const double d3 =  2.445134137142996e+00;
    static const double d4 =  3.754408661907416e+00;

    const double p_low  = 0.02425;
    const double p_high = 1.0 - p_low;

    double q, r;

    if (p < p_low) {
        q = sqrt(-2.0 * log(p));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    else if (p <= p_high) {
        q = p - 0.5;
        r = q * q;
        return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
               (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
    }
    else {
        q = sqrt(-2.0 * log(1.0 - p));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
                ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
}
// Hill's approximation for Student's t-distribution percentile
double student_t_quantile(double p, double v) {
    if (p <= 0.0 || p >= 1.0) {
        throw std::invalid_argument("Percentile must be between 0 and 1");
    }
    if (v <= 0.0) {
        throw std::invalid_argument("Degrees of freedom must be positive");
    }

    // Get the quantile for the standard normal distribution
    double z = normal_quantile(p);
    
    // Use Hill's approximation formula
    return z / sqrt((1.0 / v) + (z * z / (2.0 * v)));
}

int maxqueuelength;
void add_Qmax(int i){
    maxqueuelength = max(maxqueuelength, i);
}

double weibull_rv(double lambda, double beta){
    double u = rand()/(double)RAND_MAX; // rand() produces a uniform sample from 0 to RAND_MAX
    double w = pow(-log(u),1/beta);
    return w;
}

double pi(){ 
// well-known function in the C++ domain to get pi up to certain precision
    return atan(1)*4;
}

double trunc_normal(){ // Box - Muller Transformation
    double v = rand()/(double)RAND_MAX;
    double u = rand()/(double)RAND_MAX;
    double z = sqrt(-2*log(u))*cos(2*pi()*v);
    // We truncate the normal such that values greater than 1.5 are capped to 1.5
    // and values smaller that 0.5 are set to 0.5.
    double mean = 1.0;
    double stddev = sqrt(0.1); // Standard deviation is sqrt(0.1)
    z = mean + stddev * z;
    if(z > 1.5)
        z = 1.5;
    else if(z < 0.5) 
        z = 0.5;
    return z;
}


void run_days_regenerative_method(){
    // REGENERATIVE METHOD

    double simulation_clock = 0.;
    int ARRIVAL = 1; // arrival of a boat
    int DEPARTURE = 2; // departure of a boat
    // we also use DEPARTURE + 1,
    // which signals the case that a boat that was using two cranes departs

    queue<double> fifoqueue;
    double totalNumberOfArrivals = 0.;
    double totalWaitingTimeOfCustomers = 0.;
    double totalUnloadingTime = 0.;
    double maximumWaitingTimeOfCustomers = 0.;
    double totalWaitingTimeOfCustomersPlusUnloadingTime = 0.;
    double maximumWaitingTimeOfCustomersPlusUnloadingTime = 0.;
    int maximumQueueLength = 0;
    
    int R = int(1e6);
    vector<double> I(R); // average queue lengths for each cycle
    vector<double> T(R); // cycle lengths


    int currentNumberOfCranesBusy = 0;
    double beta = 0.709426714391227;
    double firstArrivalTime = weibull_rv(1,beta);
    PQ events;
    events.push(make_pair(firstArrivalTime, ARRIVAL));
    int c = 2;

    int cycles = 0; // currently we are at cycle 0
    int lenQ = 0; // length of the queue
    while(cycles < R){ // 90 days
        //cout << events.top().first << " " << events.top().second << endl;
        // handle next event
        double nextEventTime = events.top().first;
        int nextEventType = events.top().second;
        events.pop();
        simulation_clock = nextEventTime;

        lenQ += fifoqueue.size();

        if(currentNumberOfCranesBusy == 0 & fifoqueue.size() == 0){
            // conditions such that there are no boats at the harbor
            if(cycles == 0){
                I[cycles] = lenQ;
                T[cycles] = simulation_clock;
                cycles += 1;
            }
            else{
                I[cycles] = lenQ - I[cycles-1];
                T[cycles] = simulation_clock - T[cycles-1];
                cycles += 1;
            }
        }
        

        if(nextEventType == ARRIVAL){
            // handle arrival
            totalNumberOfArrivals += 1;
            if(currentNumberOfCranesBusy < c){
                if(currentNumberOfCranesBusy == 0 & fifoqueue.size() == 0){ 
                // this is the condition that there are no boats at the harbor
                    currentNumberOfCranesBusy = 2; // the boat uses two cranes
                    double serviceTime = trunc_normal()/2.;
                    events.push(make_pair(simulation_clock+serviceTime, DEPARTURE+1)); 
                    // departure + 1 means that the boat that leaves was using two cranes
                    totalUnloadingTime += serviceTime;
                    totalWaitingTimeOfCustomersPlusUnloadingTime += serviceTime;
                    maximumWaitingTimeOfCustomersPlusUnloadingTime = max(maximumWaitingTimeOfCustomersPlusUnloadingTime, serviceTime);
                }
                else if(fifoqueue.size() == 0){ 
                // we handle the queues in the departure of the boats
                // in this case we handle the situation in which there arrives a new boat
                // and there is a free crane
                    double serviceTime = trunc_normal();
                    currentNumberOfCranesBusy += 1;
                    events.push(make_pair(simulation_clock+serviceTime, DEPARTURE));
                    totalUnloadingTime += serviceTime;
                    totalWaitingTimeOfCustomersPlusUnloadingTime += serviceTime;
                    maximumWaitingTimeOfCustomersPlusUnloadingTime = max(maximumWaitingTimeOfCustomersPlusUnloadingTime, serviceTime);
                }
                // notice that the case that currentNumberOfCranesBusy < c 
                // and fifoqueue.size() > 0 does not happen by construction
            }
            else{ // currentNumberOfCranesBusy == c
                fifoqueue.push(nextEventTime); // we add it to the FIFO Queue
            }
            double nextArrivalTime = simulation_clock + weibull_rv(1,beta);
            events.push(make_pair(nextArrivalTime, ARRIVAL));
        }
        else{
            if(nextEventType == DEPARTURE+1){
                    currentNumberOfCranesBusy -=2;
            }
            else{
                currentNumberOfCranesBusy -=1;
            }
            // handle departure
            if(fifoqueue.size() > 0){
                // there is a customer (or more) queued                
                // we need to account for the situation in which the boat uses two cranes
                // hence, we fill the cranes given by the customer order in the FIFO queue
                int n = fifoqueue.size();
                maximumQueueLength = max(maximumQueueLength, n);

                while(fifoqueue.size() > 0 & currentNumberOfCranesBusy < c){
                    // put into service and schedule departure
                    /*
                    cout << "INFORMATION ABOUT FILLING OF THE QUEUE:" << endl;
                    cout << "currentNumberOfCranesBusy: " << currentNumberOfCranesBusy << endl;
                    cout << "fifoqueue.size(): " << fifoqueue.size() << endl;
                    cout << "fifoqueue.front(): " << fifoqueue.front() << endl;
                    */
                    currentNumberOfCranesBusy += 1;
                    double nextCustomer = fifoqueue.front();
                    fifoqueue.pop();
                    
                    totalWaitingTimeOfCustomers += (simulation_clock - nextCustomer);
                    
                    maximumWaitingTimeOfCustomers = max(maximumWaitingTimeOfCustomers,
                                                    simulation_clock - nextCustomer); 
                    
                    double nextServiceTime = trunc_normal();
                    
                    totalWaitingTimeOfCustomersPlusUnloadingTime += 
                    (simulation_clock - nextCustomer + nextServiceTime);
                    
                    maximumWaitingTimeOfCustomersPlusUnloadingTime =
                    max(maximumWaitingTimeOfCustomersPlusUnloadingTime, 
                        simulation_clock - nextCustomer + nextServiceTime);
                    
                    events.push(make_pair(simulation_clock+nextServiceTime, DEPARTURE));
                }
            }
        }
        /*
        cout << "Day: " << simulation_clock << endl;
        
        cout << "Total number of arrivals: " 
        << totalNumberOfArrivals << endl;
        
        cout << "Total waiting time of customers: " 
        << totalWaitingTimeOfCustomers << endl;
        
        cout << "Total unloading time: " 
        << totalUnloadingTime << endl;
        
        cout << "Current number of cranes busy: " 
        << currentNumberOfCranesBusy << endl;
        
        cout << "Current queue length: " 
        << fifoqueue.size() << endl;
        
        cout << "Maximum waiting time of customers: " 
        << maximumWaitingTimeOfCustomers << endl;
        
        cout << "Maximum waiting time of customers plus unloading time: " 
        << maximumWaitingTimeOfCustomersPlusUnloadingTime << endl;
        
        cout << "Total waiting time of customers plus unloading time: " 
        << totalWaitingTimeOfCustomersPlusUnloadingTime << endl;
        
        cout << "Maximum queue length: " << 
        maximumQueueLength << endl;*/
    }
    // in the same notation that in the report,
    // s_11, s_22, s_12 are the sample covariances, and
    // eta_R^2  = 1/\bar{T} (s_11 + (\bar{I}/\bar{T})^2 s_22 - 2 \bar{I}/\bar{T} s_12)
    double barT = 0.;
    double barI = 0.;
    for(int i = 0; i < R; i++){
        barT += T[i];
        barI += I[i];
    }
    barT = barT/R;
    barI = barI/R;
    double s_11 = 0.;
    double s_22 = 0.;
    double s_12 = 0.;
    for(int i = 0; i < R; i++){
        s_11 += (I[i] - barI)*(I[i] - barI);
        s_12 += (T[i] - barT)*(I[i] - barI);
        s_22 += (T[i] - barT)*(T[i] - barT);
    }
    s_11 = s_11/(R-1);
    s_22 = s_22/(R-1);
    s_12 = s_12/(R-1);
    cout << "Sample covariances:" << endl;
    cout << "s_11: " << s_11 << endl;
    cout << "s_22: " << s_22 << endl;
    cout << "s_12: " << s_12 << endl;
    cout << "Mean of T: " << barT << endl;
    cout << "Mean of I: " << barI << endl;
    double barr = barI/(double)barT;
    cout << "Mean of R (=I/T) (queue length): " << barr << endl;
    double eta = (s_11 + barr*barr*s_22 - 2*barr*s_12)/(double)barT;
    double negCI_95 = barr - 1.96*sqrt(eta)/sqrt(R);
    double posCI_95 = barr + 1.96*sqrt(eta)/sqrt(R);

    cout << "95 confidence interval: [" << negCI_95 << ", " << posCI_95 << "]" << endl;

    
    return;
}

void run_days_batch_means_method(){
    // Method of Batch Means

    double simulation_clock = 0.;
    int ARRIVAL = 1; // arrival of a boat
    int DEPARTURE = 2; // departure of a boat
    // we also use DEPARTURE + 1,
    // which signals the case that a boat that was using two cranes departs
    int t = 1e6;
    queue<double> fifoqueue;
    double totalNumberOfArrivals = 0.;
    double totalWaitingTimeOfCustomers = 0.;
    double totalUnloadingTime = 0.;
    double maximumWaitingTimeOfCustomers = 0.;
    double totalWaitingTimeOfCustomersPlusUnloadingTime = 0.;
    double maximumWaitingTimeOfCustomersPlusUnloadingTime = 0.;
    int maximumQueueLength = 0;
    
    int R = int(25); // between 5 and 30 as recommended by AG.
    vector<double> T(R);
    vector<double> Y_r(R);


    int currentNumberOfCranesBusy = 0;
    double beta = 0.709426714391227;
    double firstArrivalTime = weibull_rv(1,beta);
    PQ events;
    events.push(make_pair(firstArrivalTime, ARRIVAL));
    int c = 2;

    int lenQ = 0; // length of the queue
    while(events.top().first < t){ // t days
        //cout << events.top().first << " " << events.top().second << endl;
        // handle next event
        double nextEventTime = events.top().first;
        int nextEventType = events.top().second;
        events.pop();
        simulation_clock = nextEventTime;
        int r = int(simulation_clock / (t / R));
        Y_r[r] += fifoqueue.size();

        lenQ += fifoqueue.size();
        if(nextEventType == ARRIVAL){
            // handle arrival
            totalNumberOfArrivals += 1;
            if(currentNumberOfCranesBusy < c){
                if(currentNumberOfCranesBusy == 0 & fifoqueue.size() == 0){ 
                // this is the condition that there are no boats at the harbor
                    currentNumberOfCranesBusy = 2; // the boat uses two cranes
                    double serviceTime = trunc_normal()/2.;
                    events.push(make_pair(simulation_clock+serviceTime, DEPARTURE+1)); 
                    // departure + 1 means that the boat that leaves was using two cranes
                    totalUnloadingTime += serviceTime;
                    totalWaitingTimeOfCustomersPlusUnloadingTime += serviceTime;
                    maximumWaitingTimeOfCustomersPlusUnloadingTime = max(maximumWaitingTimeOfCustomersPlusUnloadingTime, serviceTime);
                }
                else if(fifoqueue.size() == 0){ 
                // we handle the queues in the departure of the boats
                // in this case we handle the situation in which there arrives a new boat
                // and there is a free crane
                    double serviceTime = trunc_normal();
                    currentNumberOfCranesBusy += 1;
                    events.push(make_pair(simulation_clock+serviceTime, DEPARTURE));
                    totalUnloadingTime += serviceTime;
                    totalWaitingTimeOfCustomersPlusUnloadingTime += serviceTime;
                    maximumWaitingTimeOfCustomersPlusUnloadingTime = max(maximumWaitingTimeOfCustomersPlusUnloadingTime, serviceTime);
                }
                // notice that the case that currentNumberOfCranesBusy < c 
                // and fifoqueue.size() > 0 does not happen by construction
            }
            else{ // currentNumberOfCranesBusy == c
                fifoqueue.push(nextEventTime); // we add it to the FIFO Queue
            }
            double nextArrivalTime = simulation_clock + weibull_rv(1,beta);
            events.push(make_pair(nextArrivalTime, ARRIVAL));
        }
        else{
            if(nextEventType == DEPARTURE+1){
                    currentNumberOfCranesBusy -=2;
            }
            else{
                currentNumberOfCranesBusy -=1;
            }
            // handle departure
            if(fifoqueue.size() > 0){
                // there is a customer (or more) queued                
                // we need to account for the situation in which the boat uses two cranes
                // hence, we fill the cranes given by the customer order in the FIFO queue
                int n = fifoqueue.size();
                maximumQueueLength = max(maximumQueueLength, n);

                while(fifoqueue.size() > 0 & currentNumberOfCranesBusy < c){
                    // put into service and schedule departure
                    /*
                    cout << "INFORMATION ABOUT FILLING OF THE QUEUE:" << endl;
                    cout << "currentNumberOfCranesBusy: " << currentNumberOfCranesBusy << endl;
                    cout << "fifoqueue.size(): " << fifoqueue.size() << endl;
                    cout << "fifoqueue.front(): " << fifoqueue.front() << endl;
                    */
                    currentNumberOfCranesBusy += 1;
                    double nextCustomer = fifoqueue.front();
                    fifoqueue.pop();
                    
                    totalWaitingTimeOfCustomers += (simulation_clock - nextCustomer);
                    
                    maximumWaitingTimeOfCustomers = max(maximumWaitingTimeOfCustomers,
                                                    simulation_clock - nextCustomer); 
                    
                    double nextServiceTime = trunc_normal();
                    
                    totalWaitingTimeOfCustomersPlusUnloadingTime += 
                    (simulation_clock - nextCustomer + nextServiceTime);
                    
                    maximumWaitingTimeOfCustomersPlusUnloadingTime =
                    max(maximumWaitingTimeOfCustomersPlusUnloadingTime, 
                        simulation_clock - nextCustomer + nextServiceTime);
                    
                    events.push(make_pair(simulation_clock+nextServiceTime, DEPARTURE));
                }
            }
        }
    }
    double barY_r = 0.;
    for(int i = 0; i < R; i++){
        barY_r += Y_r[i]/(t/R);
    }
    barY_r = barY_r/R;
    
    double s = 0.;
    for(int i = 0; i < R; i++){
        s += (Y_r[i]/(t/R) - barY_r)*(Y_r[i]/(t/R) - barY_r);
    }
    s = s/(R-1);
    cout << "Sample covariance: " << s << endl;
    cout << "Mean of Y (queue length): " << barY_r << endl;
    double t_value = 2.0639;//student_t_quantile(0.975, R-1 = 25);
    double negCI_95 = barY_r - t_value*sqrt(s)/sqrt(R);
    double posCI_95 = barY_r + t_value*sqrt(s)/sqrt(R);
    cout << "95 confidence interval: [" << negCI_95 << ", " << posCI_95 << "]" << endl;

    
    return;
}


int main(){
    // we also time our code
    clock_t start = clock();
    // REGENERATIVE METHOD
    // we will count the regenerative cycles,
    // and we will count R regenerative cycles,
    // instead of just doing N simulations
    unsigned int seed = 1234;
    srand(seed);
    run_days_regenerative_method();
    
    // BATCH MEANS METHOD
    // we let t -> \infty instead of t = 90 and
    // we split the space of times in R
    srand(seed);
    run_days_batch_means_method();
    
    clock_t end = clock();
    cout << "Time taken: " << ((double)(end - start))/CLOCKS_PER_SEC << " seconds" << endl;
}