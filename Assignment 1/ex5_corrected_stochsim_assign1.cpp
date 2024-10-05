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

pair<pair<double,double>, pair<double,double>> days_run_90(){
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
    int currentNumberOfCranesBusy = 0;
    double beta = 0.709426714391227;
    double firstArrivalTime = weibull_rv(1,beta);
    PQ events;
    events.push(make_pair(firstArrivalTime, ARRIVAL));
    int c = 2;
    while(events.top().first < 90){ // 90 days
        //cout << events.top().first << " " << events.top().second << endl;
        // handle next event
        double nextEventTime = events.top().first;
        int nextEventType = events.top().second;
        events.pop();
        simulation_clock = nextEventTime;

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
                    
                    cout << "INFORMATION ABOUT FILLING OF THE QUEUE:" << endl;
                    cout << "currentNumberOfCranesBusy: " << currentNumberOfCranesBusy << endl;
                    cout << "fifoqueue.size(): " << fifoqueue.size() << endl;
                    cout << "fifoqueue.front(): " << fifoqueue.front() << endl;
                    
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
        
        cout << "Day: " << simulation_clock << endl;
        
        /*cout << "Total number of arrivals: " 
        << totalNumberOfArrivals << endl;
        
        cout << "Total waiting time of customers: " 
        << totalWaitingTimeOfCustomers << endl;
        
        cout << "Total unloading time: " 
        << totalUnloadingTime << endl;
        */
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
        
        /*cout << "Maximum queue length: " << 
        maximumQueueLength << endl;
    */}
    add_Qmax(maximumQueueLength);
    return make_pair(
           make_pair(
            maximumWaitingTimeOfCustomers, 
            totalWaitingTimeOfCustomers/totalNumberOfArrivals)
           ,
           make_pair(
            maximumWaitingTimeOfCustomersPlusUnloadingTime,
            totalWaitingTimeOfCustomersPlusUnloadingTime/totalNumberOfArrivals)
           );
}

int main(){
    // we also time our code
    clock_t start = clock();
    int num_simulations = 1;
    vector<double> maxs(num_simulations);
    vector<double> means(num_simulations);
    vector<double> maxs_pls_unloading(num_simulations);
    vector<double> means_pls_unloading(num_simulations);
    pair<pair<double,double>, pair<double,double>> result;
    unsigned int seed = 1234;
    srand(seed);
    for(int i =0; i < num_simulations; ++i){
        result = days_run_90();
        maxs[i] = result.first.first;
        means[i] = result.first.second;
        maxs_pls_unloading[i] = result.second.first;
        means_pls_unloading[i] = result.second.second;
        // cout << maxs[i] << " " << means[i] << endl;
    }
    // we now compute the 95% confidence intervals for the maximum and the mean
    double alpha = 0.05;
    cout << "Maximum queue length is: " << maxqueuelength << endl;
    sort(maxs.begin(), maxs.end());
    sort(means.begin(), means.end());
    sort(maxs_pls_unloading.begin(), maxs_pls_unloading.end());
    sort(means_pls_unloading.begin(), means_pls_unloading.end());
    
    int q025 = int(num_simulations*alpha/2);
    int q975 = int(num_simulations*(1-alpha/2));
    
    double max_CI = maxs[q025], max_CI2 = maxs[q975];
    double mean_CI = means[q025], mean_CI2 = means[q975];
    double max_CI_pls_unloading = 
    maxs_pls_unloading[q025], max_CI2_pls_unloading = maxs_pls_unloading[q975];
    double mean_CI_pls_unloading = 
    means_pls_unloading[q025], mean_CI2_pls_unloading = means_pls_unloading[q975];

    cout << "The maximum waiting time of customers is between " 
    << max_CI << " and " << max_CI2 << endl;
    cout << "The mean waiting time of customers is between " 
    << mean_CI << " and " << mean_CI2 << endl;
    cout << "The maximum waiting time understood as waiting time until unloaded is between "
    << max_CI_pls_unloading << " and " << max_CI2_pls_unloading << endl;
    cout << "The mean waiting time understood as waiting time until unloaded is between "
    << mean_CI_pls_unloading << " and " << mean_CI2_pls_unloading << endl;
    clock_t end = clock();
    cout << "Time taken: " 
    << ((double)(end - start))/CLOCKS_PER_SEC << " seconds" << endl;
}