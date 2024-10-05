from numpy import random
import math
import matplotlib.pyplot as plt
from scipy.stats import invgauss

###Exercise 5:

def normal_generator(mu, sigma):
    return float(random.normal(mu, sigma, 1)[0]) ###to get a number and not an array

def weibull_generator(beta):
    return float(random.weibull(beta, 1)[0]) ###to get a number and not an array

def add_event(eventList, event,t):
    eventList[t] = event

def getNextEvent(eventList):
    firstKey = sorted(eventList.keys())[0]  
    event = eventList.pop(firstKey)  
    return [firstKey, event]

def nextBoatFromQueue(queue):
    firstKey = sorted(queue.keys())[0]  
    boat = queue.pop(firstKey)
    return [firstKey, boat]

def arrivalEvent(w, n, l,t, eventList, queue, mu, sigma2, beta, o):
    if n == 0 and len(queue) == 0 and o == 0:
        unloading = normal_generator(mu, sigma2)
        if unloading < 0.5: unloading = 0.5
        elif unloading > 1.5: unloading = 1.5
        unloading /= 2
        nextDeparture = unloading + t
        n = 1
        o = 2
        w += unloading
        add_event(eventList, 'DEPARTURE' ,nextDeparture)
    elif n == 1 and len(queue) == 0 and o == 1:
        unloading = normal_generator(mu, sigma2)
        if unloading < 0.5: unloading = 0.5
        elif unloading > 1.5: unloading = 1.5
        nextDeparture = unloading + t
        n = 2
        o = 2
        w += unloading
        add_event(eventList, 'DEPARTURE' ,nextDeparture)
    else:
        queue[t] = l ###new boat in the line
    nextArrival = t + weibull_generator(beta)
    l += 1
    add_event(eventList, "ARRIVAL", nextArrival)
    return w, n, l, o, queue

def departureEvent(w, n, mu, sigma2, queue, t, maxW, eventList, o):
    if len(queue) == 1 and n == 1 and o >= 1:
        arrivalTimeOfBoat = nextBoatFromQueue(queue)[0]
        unloading = normal_generator(mu, sigma2)
        if unloading < 0.5: unloading = 0.5
        elif unloading > 1.5: unloading = 1.5
        unloading /= 2
        totalWaitingTime = t - arrivalTimeOfBoat + unloading 
        w += totalWaitingTime
        waitingOnService = t - arrivalTimeOfBoat
        if waitingOnService > maxW: maxW = waitingOnService ##might do a list of waiting times and return the max
        nextDeparture = t + unloading
        add_event(eventList, "DEPARTURE", nextDeparture)
        n = 1
        o = 2
    elif len(queue) >= 1 and n == 2 and o == 2:
        arrivalTimeOfBoat = nextBoatFromQueue(queue)[0]
        unloading = normal_generator(mu, sigma2)
        if unloading < 0.5: unloading = 0.5
        elif unloading > 1.5: unloading = 1.5
        totalWaitingTime = t - arrivalTimeOfBoat + unloading 
        w += totalWaitingTime
        waitingOnService = t - arrivalTimeOfBoat
        if waitingOnService > maxW: maxW = waitingOnService ##might do a list of waiting times and return the max
        nextDeparture = unloading + t
        add_event(eventList, "DEPARTURE", nextDeparture)
        n = 2
        o = 2
    elif len(queue) > 1 and n == 1 and o > 0 :
        arrivalTimeOfBoat1 = nextBoatFromQueue(queue)[0]
        arrivalTimeOfBoat2 = nextBoatFromQueue(queue)[0]
        unloading1 = normal_generator(mu, sigma2)
        unloading2 = normal_generator(mu, sigma2)
        if unloading1 < 0.5: unloading1 = 0.5
        if unloading2 < 0.5: unloading2 = 0.5
        if unloading1 > 1.5: unloading1 = 1.5
        if unloading2 > 1.5: unloading2 = 1.5
        totalWaitingTime1 = t - arrivalTimeOfBoat1 + unloading1 ##check
        totalWaitingTime2 = t - arrivalTimeOfBoat2 + unloading2 ##check
        w += totalWaitingTime1 + totalWaitingTime2
        waitingOnService1 = t - arrivalTimeOfBoat1
        waitingOnService2 = t - arrivalTimeOfBoat2
        if max(waitingOnService1, waitingOnService2) > maxW: maxW = max(waitingOnService1, waitingOnService2)
        nextDeparture1 = t + unloading1
        nextDeparture2 = t + unloading2
        add_event(eventList, "DEPARTURE", nextDeparture1)
        add_event(eventList, "DEPARTURE", nextDeparture2)
        n = 2
        o = 2
    elif len(queue) == 0 and o == 2 and n == 1:
        n = 0
        o = 0
    elif len(queue) == 0 and o == 1:
        n = 0
        o = 0
    elif len(queue) == 0 and o == 2 and n == 2:
        n = 1
        o = 1
    return w, maxW, n, o, queue

def nextEvent(eventList, t, w, n, l, queue, mu, sigma2, beta, o, maxW, q):
    event = getNextEvent(eventList)
    eventTime = event[0]
    eventType = event[1]
    t = eventTime
    if eventType == "ARRIVAL":
        w, n, l, o, queue = arrivalEvent(w, n, l, t, eventList, queue, mu, sigma2, beta, o)
    elif eventType == "DEPARTURE":
        w, maxW, n, o, queue = departureEvent(w, n, mu, sigma2, queue, t, maxW, eventList, o)
    q = max(q, len(queue))
    return [t, w, maxW, n, l, o, eventType,q]

def simulation_discrete_events():
    t = 0 ##simulation clock

    c = 2 ##number of cranes
    b = 2 ##number of berths

    beta = 0.7094267 ##calculated with WolframAlpha
    mu = 1
    sigma2 = 0.1
    
    l = 0 ##number of boats that came in the harbor
    n = 0 ##number of occupied berths
    w = 0 ##waiting time
    o = 0 ##number of occupied cranes
    maxW = 0 ##maximum waiting time

    numberOfDays = 90 ###number of days the harbor is operating
    eventList = {} ##type time: event type (departure, arrival)
    queue = {} ##type time of arrival: i 
    q = 0 ###max length of queue (for test)

    firstArrival = weibull_generator(beta)
    add_event(eventList, "ARRIVAL", firstArrival)
    t = firstArrival
    while t <= numberOfDays:
        [t, w, maxW, n, l, o, eventType, q]= nextEvent(eventList, t, w, n, l, queue, mu, sigma2, beta, o, maxW, q)
        ##print(eventType + " at time " + str(t))
        ##print("len of queue: " + str(len(queue)))
        ##print("occupied berths: " + str(n) + " occupied cranes: " + str(o))

    averageW = w / l

    print("The total waiting time: " + str(w))
    print("The average waiting time: " + str(averageW))
    print("Max waiting time: " + str(maxW))

    return [w, averageW , maxW, t, l, q]

""" beta = 0.7094267
##print(simulation_discrete_events())

eventList = {1.4:"ARRIVAL", 5.2: "DEPARTURE", 0.4: "DEPARTURE"}
queue = {7.3: 1, 7.7: 2, 9.3: 3}

w = 2
maxW = 4
w, maxW, n, o = departureEvent(w, 2, 1, 0.1, queue, 7, maxW, eventList,2)
##getNextEvent(eventList)
print(eventList)
print(w)
print(maxW)
print(n)
print(o)  """

def repetitions(n = 10000):
    N = 0
    L = 0
    M = 0
    W = 0
    Q = 0
    w = 0
    averageW = 0
    j = 0
    maxW = 0
    for i in range(n):
        w, averageW, maxW, t, l, q = simulation_discrete_events()
        N += l
        W += w
        M = max(M, maxW)
        L += l
        if max(q, Q) == q:
            w = w
            Q = q
            averageW = averageW
            j = i
            maxW = maxW
    average = W/ L
    print("At day " + str(j) + " we have total waiting time: " + str(w) + " , average waiting time: " + str(averageW) + ", max waiting time: " + str(maxW) + " and max length of queue: " + str(q)  )
    return (W, average, M, Q)

##simulation_discrete_events()
print(repetitions())