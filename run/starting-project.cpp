// Copyright © 2021 Giorgio Audrito. All Rights Reserved.


// [INTRODUCTION]
//! Importing the FCPP library.
#include "lib/fcpp.hpp"

const int NODE_NUM = 300;

namespace fcpp {

namespace coordination {

namespace tags {
    struct node_color {};
    struct node_size {};
    struct node_shape {};

    struct flow_field {};
    struct obstruction {};

    //fields of outgoing capacities
    struct capacity_field {};

    struct node_distance_to_sink{};

    struct out_flow{};
    struct in_flow{};

    struct obstruction_condition{};


}

//! @brief The maximum communication range between nodes.
constexpr size_t communication_range = 100;

const size_t dim = 2;



/*
In the description of the algorithm, weighted directed graphs on the set of devices will be 
identified with fields of reals:
if  d1 is a device, d2 a neighbour of d1 and the field F in d1  has value w nonzero in d2,
then we say that the graph associated with F has an arc from d1 to d2 with weight w.

This algorithm is inspired to the Ford–Fulkerson method: as long as there is an admissible 
path from source to sink in the residual graph, we push flow along the path.
We detect the presence of such paths simply by considering the distances from source and 
to sink along the residual graph.
If a node has both finite distance from source and finite distance to sink, then there is an
admissible path passing through it.
*/


/* Here we implement our rules to assign dynamically the capacities between nodes.
These capacities will remain unaffected by the updates in our flow.
*/
FUN field<real_t> capacity_v0(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid!=id? 1.0 : 0.0 ;}, ids);
}
/*
FUN field<real_t> capacity_v1(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ;}, ids);
}
*/


FUN field<real_t> capacity_v2(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ? id-node.uid:node.uid-id;}, ids);
}

/*
FUN field<real_t> capacity_v3(ARGS){ CODE
    field<device_t> ids = nbr_uid(CALL);
    return map_hood([&](device_t id){ return node.uid<id ? id-node.uid:0;}, ids);
}
*/


// Rough method to switch between capacities
FUN field<real_t> capacity(ARGS){ CODE
    
    return capacity_v0(CALL);
}


// Usual definition for the residual graph. Notice that it has only nonnegative weights.
// As already mentioned, it can have cycles.
FUN field<real_t> residual_capacity(ARGS, field<real_t> flow){ CODE
    //return mux(flow>0, 0.0, capacity(CALL)+flow );
    return capacity(CALL)-flow ;
}

//This is just another function to sum up the values in a field
real_t sum(field<real_t> input){

    real_t tmp= 0.0;
    for (details::field_iterator<field<real_t>> it(input); !it.end(); ++it){
        tmp += it.value();
    }
    return tmp;
}

field<real_t> truncate(field<real_t> inputField, real_t inputValue){

    return inputValue >= 0
    ? map_hood([&](real_t i){
                                        real_t tmp= std::max(std::min(i, inputValue),0.0);
                                        inputValue-= tmp;
                                        return tmp;
                                    }, inputField)
    : -truncate(-inputField, - inputValue);
}

/*
Function that monitors the ability of a node to push or receive flow.
Strictly negative excess means that a node receives more flow in ingress so it can push some.
Source is always able to push flow.
Strictly positive excess means that a node pushes more flow so it has to receive some.
Sink is always ready to receive flow.
*/
FUN real_t excess(ARGS, field<real_t> flow){ CODE

    bool is_source = node.uid==0;
    bool is_sink = node.uid == NODE_NUM-1;
    return is_source
                ?INF
                :is_sink
                    ? -INF
                    :sum( flow);
}


// Returns the distance to the closest sink-like node
FUN real_t to_sink(ARGS, field<real_t> flow){ CODE
    bool is_sink_ = (node.uid == NODE_NUM-1);
    return abf_distance(CALL, is_sink_, [&](){return mux(capacity(CALL) + flow>0 && flow<= 0, 1.0, INF);});
}


//Updates the flow adding the increment
FUN field<real_t> update_flow(ARGS, field<real_t>& flow_){ CODE

        field<real_t> capacity_n = capacity(CALL);
        
        //safety conditions
        mod_other(CALL, flow_) = 0.0;
        field<real_t> flow = mux( flow_>0 , std::min(flow_, capacity_n), std::max(flow_, -capacity_n));
        //

        real_t excess_n = excess(CALL, flow);

        real_t to_sink_n = to_sink(CALL, flow);

        field<real_t> forward = truncate( (nbr(CALL, to_sink_n)<to_sink_n) * (capacity_n + flow),
                                     excess_n);

        return  -flow 
        + forward //push forward
        + truncate(flow, excess(CALL, flow-forward) ); //then push backward
}





//! @brief Main function.
MAIN() {
    // import tag names in the local scope.
    using namespace tags;

    // References
    field<real_t>& capacity_ = node.storage(capacity_field{});
    field<real_t>& flow_ = node.storage(flow_field{});
    real_t& to_sink_ = node.storage(node_distance_to_sink{});
    real_t& obstruction_ = node.storage(obstruction{});
    real_t& out_flow_ = node.storage(out_flow{});
    real_t& in_flow_ = node.storage(in_flow{});
    real_t& obstruction_condition_ = node.storage(obstruction_condition{});


    // Usage of node storage
    bool is_source = node.uid==0;
    bool is_sink = node.uid == NODE_NUM-1;

    node.storage(node_shape{}) = is_source
                                    ?shape::star 
                                    :is_sink
                                        ? shape::tetrahedron
                                        : shape::sphere;

    // This is the only structure that node requires to manage
    flow_ = nbr(CALL, field<real_t>(0.0),[&](field<real_t> flow){
            return update_flow(CALL, flow);
            });
    
    node.storage(node_size{}) = 8;

    // These other structures are just aimed at monitoring the behaviour of system
    capacity_ = capacity(CALL);
    
    to_sink_ = to_sink(CALL, flow_);


    /*
    In this structurre we monitor how much flow source pushes and
    how much flow sink receives. Hopefully they're equal in absolute module
    */
    out_flow_= sum(mux(flow_>0, flow_, 0.0));
    in_flow_= sum(mux(flow_<0, -flow_, 0.0));


    obstruction_=  to_sink_<INF
            ? sum(capacity(CALL) + flow_)
            : -sum(capacity(CALL) + flow_);


    //eventually true
    obstruction_condition_ = obstruction_ <= old(CALL, obstruction_);



    
    
    /*
    Nodes that have flow to push are GREEN; those that need to receive flow are RED;
    Nodes that guess to be part of an admissible path from a source-like to a sink-like
    are YELLOW; other nodes are WHITE.
    */

    node.storage(node_color{}) =   sum(flow_)>0
                                    ? color(GREEN)
                                    : sum(flow_)<0
                                        ? color(RED)
                                        : node.storage(out_flow{})>0
                                            ?color(YELLOW)
                                            : color(BLACK);
}
//! @brief Export types used by the main function.
FUN_EXPORT main_t = export_list<device_t, field<real_t>, real_t, bool, field<int>>;

} // namespace coordination

// [SYSTEM SETUP]

//! @brief Namespace for component options.
namespace option {

//! @brief Import tags to be used for component options.
using namespace component::tags;
//! @brief Import tags used by aggregate functions.
using namespace coordination::tags;

//! @brief Number of people in the area.
constexpr int node_num = NODE_NUM;
//! @brief Dimensionality of the space.
constexpr size_t dim = 2;

//! @brief Description of the round schedule.
using round_s = sequence::periodic<
    distribution::interval_n<times_t, 0, 1>,    // uniform time in the [0,1] interval for start
    distribution::weibull_n<times_t, 10, 1, 10> // weibull-distributed time for interval (10/10=1 mean, 1/10=0.1 deviation)
>;
//! @brief The sequence of network snapshots (one every simulated second).
using log_s = sequence::periodic_n<1, 0, 1>;
//! @brief The sequence of node generation events (node_num devices all generated at time 0).
using spawn_s = sequence::multiple_n<node_num, 0>;
//! @brief The distribution of initial node positions.
using rectangle_d = distribution::rect_n<1, 0, 0, 350, 350>;
//! @brief The contents of the node storage as tags and associated types.
using store_t = tuple_store<
    node_color,                         color,
    node_size,                          double,
    node_shape,                         shape,
    node_distance_to_sink,              real_t,
    obstruction,                        real_t,
    capacity_field,                     field<real_t>,
    flow_field,                         field<real_t>,
    out_flow,                           real_t,
    in_flow,                            real_t,
    obstruction_condition,              real_t
>;
//! @brief The tags and corresponding aggregators to be logged (change as needed).
using aggregator_t = aggregators< 
    out_flow,                   aggregator::max<real_t>,
    in_flow,                    aggregator::max<real_t>,
    obstruction,                aggregator::sum<real_t>,
    obstruction_condition,      aggregator::min<real_t>
>;

//! @brief The general simulation options.
DECLARE_OPTIONS(list,
    parallel<true>,      // multithreading enabled on node rounds
    synchronised<false>, // optimise for asynchronous networks
    program<coordination::main>,   // program to be run (refers to MAIN above)
    exports<coordination::main_t>, // export type list (types used in messages)
    retain<metric::retain<2,1>>,   // messages are kept for 2 seconds before expiring
    round_schedule<round_s>, // the sequence generator for round events on nodes
    log_schedule<log_s>,     // the sequence generator for log events on the network
    spawn_schedule<spawn_s>, // the sequence generator of node creation events on the network
    store_t,       // the contents of the node storage
    aggregator_t,  // the tags and corresponding aggregators to be logged
    init<
        x,      rectangle_d // initialise position randomly in a rectangle for new nodes
    >,
    dimension<dim>, // dimensionality of the space
    connector<connect::fixed<100, 1, dim>>, // connection allowed within a fixed comm range
    //label_text_tag<node_text>, // the text to associate to a node
    label_size_val<1, 4>,        // the size of the text
    shape_tag<node_shape>, // the shape of a node is read from this tag in the store
    size_tag<node_size>,   // the size  of a node is read from this tag in the store
    color_tag<node_color>  // the color of a node is read from this tag in the store
);

} // namespace option

} // namespace fcpp


//! @brief The main function.
int main() {
    using namespace fcpp;

    //! @brief The network object type (interactive simulator with given options).
    using net_t = component::interactive_simulator<option::list>::net;
    //! @brief The initialisation values (simulation name).
    auto init_v = common::make_tagged_tuple<option::name>("Starting Project");
    //! @brief Construct the network object.
    net_t network{init_v};
    //! @brief Run the simulation until exit.
    network.run();
    return 0;
}
