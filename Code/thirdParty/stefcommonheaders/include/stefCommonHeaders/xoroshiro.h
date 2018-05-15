#pragma once
#include <stdint.h>
#include <assert.h>


/**************************************************************************************************
 *                          THE XOROSHIRO 128 + RANDOM NUMBER GENERATOR                           *
 **************************************************************************************************/

/*
 * Shamelessly stolen from http://xoroshiro.di.unimi.it/xoroshiro128plus.c (creative commons C0)
 * All credit is due to Sebastiano Vigna and David Blackman
 *
 */

/* This is the successor to xorshift128+. It is the fastest full-period
   generator passing BigCrush without systematic failures, but due to the
   relatively short period it is acceptable only for applications with a
   mild amount of parallelism; otherwise, use a xorshift1024* generator.

   Beside passing BigCrush, this generator passes the PractRand test suite
   up to (and included) 16TB, with the exception of binary rank tests,
   which fail due to the lowest bit being an LFSR; all other bits pass all
   tests. We suggest to use a sign test to extract a random Boolean value.

   Note that the generator uses a simulated rotate operation, which most C
   compilers will turn into a single instruction. In Java, you can use
   Long.rotateLeft(). In languages that do not make low-level rotation
   instructions accessible xorshift128+ could be faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */


const static uint64_t Xoroshiro_MinimumGeneratedNumber = 0;
const static uint64_t Xoroshiro_MaximumGeneratedNumber = UINT64_MAX;


typedef struct Xoroshiro_State_Struct {
    uint64_t lower;
    uint64_t upper;
} Xoroshiro_State_t;

static inline
uint64_t
Xoroshiro_RotateLeft(
    const uint64_t x,
    int k
) {
    return (x << k) | (x >> (64 - k));
}

// Returns a random uint64_t
static inline
uint64_t
Xoroshiro_GetNextNumber(
    Xoroshiro_State_t* state
) {
    const uint64_t s0 = state->lower;
    uint64_t s1 = state->upper;
    const uint64_t result = s0 + s1;

    s1 ^= s0;
    state->lower = Xoroshiro_RotateLeft(s0, 55) ^ s1 ^ (s1 << 14); // a, b
    state->upper = Xoroshiro_RotateLeft(s1, 36); // c

    return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

static inline
void
Xoroshiro_JumpNumbers(
    Xoroshiro_State_t* state
) {
    static const uint64_t JUMP[] = { 0x8a5cd789635d2dff, 0x121fd2155c472f96 };

    uint64_t newLowerState = 0;
    uint64_t newUpperState = 0;
    for (int i = 0; i < (int)(sizeof(JUMP) / sizeof(*JUMP)); i++)
        for (int b = 0; b < 64; b++) {
            if (JUMP[i] & 1ULL << b) {
                newLowerState ^= state->lower;
                newUpperState ^= state->upper;
            }
            Xoroshiro_GetNextNumber(state);
        }

    state->lower = newLowerState;
    state->upper = newUpperState;
}

static inline
uint64_t
Xoroshiro_SplitMix64(
    uint64_t* x
) {
    uint64_t z = (*x += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

static inline
Xoroshiro_State_t
Xoroshiro_GetInitialState(
    uint64_t seed,
    const unsigned threadNo
) {
    Xoroshiro_State_t initialState;
    initialState.lower = seed;
    initialState.upper = Xoroshiro_SplitMix64(&seed);

    // jump to my section; inefficient if threadNo is high.
    // Then we would need to distribute the initial state
    for (unsigned j = 0; j < threadNo; ++j) {
        Xoroshiro_JumpNumbers(&initialState);
    }

    return initialState;
}

static inline
double
Xoroshiro_GenerateRealInClosedInterval(
    const double lowerBound,
    const double upperBound,
    Xoroshiro_State_t* state
) {
    assert(lowerBound <= upperBound);
    const uint64_t randomInteger = Xoroshiro_GetNextNumber(state);

    const double intervalLength = upperBound - lowerBound;
    const double randomReal = lowerBound + intervalLength * (double) (randomInteger - Xoroshiro_MinimumGeneratedNumber) / (double) (Xoroshiro_MaximumGeneratedNumber - Xoroshiro_MinimumGeneratedNumber);

    assert(lowerBound <= randomReal && randomReal <= upperBound);

    return randomReal;
}

/**************************************************************************************************
 *                                     XOROSHIRO AS A GSL RNG                                     *
 **************************************************************************************************/

#ifdef XOROSHIRO_GSL
#include <gsl/gsl_rng.h>

static inline
void
Xoroshiro_GSL_Set(
    void* state,
    unsigned long seed
) {
    *((Xoroshiro_State_t*) state) = Xoroshiro_GetInitialState(seed, 0);
}

static inline
unsigned long
Xoroshiro_GSL_Get(
    void* state
) {
    return (unsigned long) Xoroshiro_GetNextNumber((Xoroshiro_State_t*) state);
}

static inline
double
Xoroshiro_GSL_GetDouble(
    void* state
) {
    return Xoroshiro_GenerateRealInClosedInterval(0.0, 1.0, (Xoroshiro_State_t*) state);
}

static const gsl_rng_type Xoroshiro_GSL_rng_type = {
    "Xoroshiro", /*name*/
    Xoroshiro_MinimumGeneratedNumber, /*min*/
    Xoroshiro_MaximumGeneratedNumber, /*max*/
    sizeof(Xoroshiro_State_t), /*size*/
    Xoroshiro_GSL_Set, /*set*/
    Xoroshiro_GSL_Get, /*get*/
    Xoroshiro_GSL_GetDouble, /*get_double*/
};

static const gsl_rng_type* Xoroshiro_GSL_rng = &Xoroshiro_GSL_rng_type;

#endif

/**************************************************************************************************
 *             ADAPTATION OF THE XOROSHIRO RNG AS C++ 11 RANDOM NUMBER ENGINE CONCEPT             *
 **************************************************************************************************/
#ifdef __cplusplus

#include <stdexcept>
// satisfies the concept UniformRandomBitGenerator
class XoroshiroRandomNumberEngine {



    public:
        // named with underscore for compatibility with other RNGs in <random>
        using result_type = uint64_t;
        const static result_type default_seed = 42;



        inline
        XoroshiroRandomNumberEngine(
        ) {
            seed(default_seed);

        }


        inline
        XoroshiroRandomNumberEngine(
            const uint64_t seedToUse,
            const unsigned threadNo = 0
        ) {
            seed(seedToUse, threadNo);
            #ifdef XOROSHIRO_GSL
            gsl_rng_ptr_M = nullptr;
            #endif
        }

        inline
        XoroshiroRandomNumberEngine(
            const XoroshiroRandomNumberEngine& other
        ): state_M(other.state_M), seedUsed_M(other.seedUsed_M) {
            #ifdef XOROSHIRO_GSL
            gsl_rng_ptr_M = nullptr;
            #endif

        }

        inline
        XoroshiroRandomNumberEngine&
        operator=(
            const XoroshiroRandomNumberEngine& other
        ) {
            state_M = other.state_M;
            seedUsed_M = other.seedUsed_M;

            return *this;
        }

        inline
        XoroshiroRandomNumberEngine&
        operator=(
            XoroshiroRandomNumberEngine&& other
        ) {
            state_M = std::move(other.state_M);
            seedUsed_M = std::move(other.seedUsed_M);
            other.maybeFreeGSLRNG();
            return *this;
        }



        inline
        XoroshiroRandomNumberEngine(
            XoroshiroRandomNumberEngine&& other
        ): state_M(std::move(other.state_M)), seedUsed_M(std::move(other.seedUsed_M)) {
            other.maybeFreeGSLRNG();
        }


        inline
        ~XoroshiroRandomNumberEngine() {
            maybeFreeGSLRNG();
        }

        inline
        void
        seed(
            uint64_t seedToUse,
            const unsigned threadNo = 0
        ) {
            seedUsed_M = seedToUse;
            state_M = Xoroshiro_GetInitialState(seedToUse, threadNo);
        }

        inline
        result_type
        getSeed(
        ) const {
            return seedUsed_M;
        }

        inline
        Xoroshiro_State_t
        getState(
        ) const {
            return state_M;
        }

        inline static constexpr
        result_type
        min(
        ) {
            return Xoroshiro_MinimumGeneratedNumber;
        }


        inline static constexpr
        result_type
        max(
        ) {
            return Xoroshiro_MaximumGeneratedNumber;
        }

        inline
        result_type
        operator()(
        ) {
            return Xoroshiro_GetNextNumber(&state_M);
        }

        #ifdef XOROSHIRO_GSL
        inline
        operator gsl_rng*() {
            maybeAllocGSLRNG();
            return gsl_rng_ptr_M;
        }
        #endif



    private:
        Xoroshiro_State_t state_M;
        result_type seedUsed_M;

        #ifdef XOROSHIRO_GSL
        gsl_rng* gsl_rng_ptr_M;
        #endif

        inline
        void
        maybeAllocGSLRNG() {
            #ifdef XOROSHIRO_GSL
            if (gsl_rng_ptr_M != nullptr) {
                return;
            }
            gsl_rng_ptr_M = gsl_rng_alloc(Xoroshiro_GSL_rng);
            if (gsl_rng_ptr_M == nullptr) {
                throw std::runtime_error("Could not allocate Xoroshiro as a GSL RNG!");
            }
            gsl_rng_ptr_M->state = static_cast<void*>(&state_M);
            assert(gsl_rng_ptr_M->state != nullptr);
            //printf("allocated\n");
            #endif
        }

        inline
        void
        maybeFreeGSLRNG() {
            #ifdef XOROSHIRO_GSL
            if (gsl_rng_ptr_M != nullptr) {
                gsl_rng_free(gsl_rng_ptr_M);
            }
            gsl_rng_ptr_M = nullptr;
            //printf("freed\n");
            #endif
        }
};

#endif
