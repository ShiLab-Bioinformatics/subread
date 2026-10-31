/* Direction compatibility of breakpoint retained flanks in subjunc fusion mode. */
#ifndef SUBJUNC_EVENT_DIRECTION_H
#define SUBJUNC_EVENT_DIRECTION_H

/* A side value of one denotes reference bases to the right of the breakpoint. */
static inline int subjunc_event_direction_allowed(int search_back, int reversed,
        int enter_large, int small_right, int large_right, int strand_jumped)
{
    int incoming = enter_large ? !!large_right : !!small_right;
    int outgoing = enter_large ? !!small_right : !!large_right;
    int next_reversed = !!reversed ^ !!strand_jumped;
    int expected_incoming = !!reversed ^ !!search_back;
    int expected_outgoing = next_reversed ^ !search_back;
    return incoming == expected_incoming && outgoing == expected_outgoing;
}

#endif
