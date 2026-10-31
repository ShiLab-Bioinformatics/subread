#include <assert.h>
#include <stdio.h>
#include "../../src/subjunc-event-direction.h"

int main(void)
{
    int back, reversed, enter_large, small_right, large_right, jumped;
    int valid = 0, invalid = 0;
    for(back = 0; back < 2; back++)
    for(reversed = 0; reversed < 2; reversed++)
    for(enter_large = 0; enter_large < 2; enter_large++)
    for(small_right = 0; small_right < 2; small_right++)
    for(large_right = 0; large_right < 2; large_right++)
    for(jumped = 0; jumped < 2; jumped++) {
        int incoming = enter_large ? large_right : small_right;
        int outgoing = enter_large ? small_right : large_right;
        int approach_increases = back ? reversed : !reversed;
        int leaving_increases = jumped ? !approach_increases : approach_increases;
        int expected = incoming != approach_increases && outgoing == leaving_increases;
        int actual = subjunc_event_direction_allowed(back, reversed, enter_large,
                                                     small_right, large_right, jumped);
        if(actual != expected) {
            fprintf(stderr, "direction mismatch: %d %d %d %d %d %d\n", back, reversed,
                    enter_large, small_right, large_right, jumped);
            return 1;
        }
        if(actual) valid++; else invalid++;
    }
    assert(valid == 16 && invalid == 48);
    /* Ordinary junction: front entering the large side and back entering the
       small side on the primary orientation are the false back-to-back paths. */
    assert(!subjunc_event_direction_allowed(0, 0, 1, 0, 1, 0));
    assert(!subjunc_event_direction_allowed(1, 0, 0, 0, 1, 0));
    assert(subjunc_event_direction_allowed(0, 0, 0, 0, 1, 0));
    assert(subjunc_event_direction_allowed(1, 0, 1, 0, 1, 0));
    puts("subjunc event-direction tests passed");
    return 0;
}
