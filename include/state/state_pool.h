#ifndef _STATE_POOL_H
#define _STATE_POOL_H

#include <list>

#include "base/types.h"
#include "state/state.h"

class StatePool {
public:
    StatePool() = default;
    StatePool(sz_t n_reserved, const State& s) { reserve(n_reserved, s); }

    void reserve(sz_t n_reserved, const State& s) {
        for (sz_t i = 0; i < n_reserved; ++i) pool.push_front(new_similar(s));
    }

    ~StatePool() {
        for (State* s : pool) delete s;
    }

    class StateGuard {
    public:
        StateGuard(StatePool& pool, State& state) : state(state), pool(pool) {}
        ~StateGuard() { pool.recycle(state); }

        State& state;
    private:
        StatePool& pool;
    };

    StateGuard allocate_similar(const State& s) {
        for (auto it = pool.begin(); it != pool.end(); ++it) {
            State* x = *it;
            if (x->total_dims() == s.total_dims()) {
                pool.erase(it);
                return StateGuard(*this, *x);
            }
        }
        
        return StateGuard(*this, *new_similar(s));
    }

private:
    std::list<State*> pool;

    State* new_similar(const State& s) {
        State* x = new State();
        x->make_similar(s);
        return x;
    }

    friend class StateGuard;
    void recycle(State& s) {
        pool.push_front(&s);
    }
};

#endif // _STATE_POOL_H
