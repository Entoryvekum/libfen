#ifndef CPU_RELAX_H
    #define CPU_RELAX_H
    #if defined(_MSC_VER)
        #include <intrin.h>
        static inline void cpu_pause() {
            #if defined(_M_ARM64) || defined(_M_ARM)
                __yield();
            #else
                _mm_pause();
            #endif
        }
    #else
        static inline void cpu_pause() {
            #if defined(__i386__) || defined(__x86_64__)
                __asm__ __volatile__("pause");
            #elif defined(__aarch64__) || defined(__arm__)
                __asm__ __volatile__("yield");
            #else
                __asm__ __volatile__("" ::: "memory");
            #endif
        }
    #endif
#endif