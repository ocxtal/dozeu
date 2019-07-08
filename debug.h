
/**
 * @file debug.h
 */
#ifndef _DEBUG_H_INCLUDED
#define _DEBUG_H_INCLUDED

/* include global header *before* we include individual dependencies */
// #include "common.h"
// #include "xprintf.h"


#if defined(DEBUG) && !defined(NDEBUG_BLOCK)
#  define DEBUG_BLOCK
#endif
#if defined(DEBUG) && !defined(NDEBUG_PRINT)
#  define DEBUG_PRINT
#endif
#if defined(DEBUG) && !defined(NDEBUG_TRAP)
#  define DEBUG_TRAP
#endif
#if defined(DEBUG) && !defined(NDEBUG_ASSERT)
#  define DEBUG_ASSERT
#endif
#if defined(DEBUG) && !defined(NDEBUG_SCAN)
#  define DEBUG_SCAN
#endif


#if defined(NDEBUG)
#  undef DEBUG_BLOCK
#  undef DEBUG_PRINT
#  undef DEBUG_TRAP
#  undef DEBUG_ASSERT
#  undef DEBUG_SCAN
#endif

/**
 * static assertion macros
 */
#define _sa_cat_intl(x, y)		x##y
#define _sa_cat(x, y)			_sa_cat_intl(x, y)
#define _static_assert(expr)	typedef char _sa_cat(_st, __LINE__)[(expr) ? 1 : -1]

/**
 * color outputs
 */
#define RED(x)			"\x1b[31m" x "\x1b[39m"
#define GREEN(x)		"\x1b[32m" x "\x1b[39m"
#define YELLOW(x)		"\x1b[33m" x "\x1b[39m"
#define BLUE(x)			"\x1b[34m" x "\x1b[39m"
#define MAGENTA(x)		"\x1b[35m" x "\x1b[39m"
#define CYAN(x)			"\x1b[36m" x "\x1b[39m"
#define WHITE(x)		"\x1b[37m" x "\x1b[39m"

/**
 * blockwise control
 */
#undef debugblock
#ifdef DEBUG_BLOCK
#define debugblock(x) { x; }
#else
#define debugblock(x) {}
#endif

/**
 * printers
 */
#undef debug
#undef debug_impl
#undef dbprintf
#undef dump

#ifdef DEBUG_PRINT
#define debug(...) { \
	debug_impl(__VA_ARGS__, ""); \
}
#define debug_impl(fmt, ...) { \
	dbprintf("[%s: %s(%d)] " fmt "%s\n", __FILE__, __func__, __LINE__, __VA_ARGS__); \
}
#define dbprintf(fmt, ...) { \
	fprintf(stderr, fmt, __VA_ARGS__); \
}
/* compatible with dump in unittest.h */
#if 0
#define dump(ptr, len) ({ \
	uint64_t size = (((len) + 15) / 16 + 1) * \
		(strlen("0x0123456789abcdef:") + 16 * strlen(" 00a") + strlen("  \n+ margin")) \
		+ strlen(#ptr) + strlen("\n`' len: 100000000"); \
	uint8_t *_ptr = (uint8_t *)(ptr); \
	char *_str = alloca(size); \
	char *_s = _str; \
	/* make header */ \
	_s += sprintf(_s, "\n`%s' len: %" PRId64 "\n", #ptr, (int64_t)len); \
	_s += sprintf(_s, "                   "); \
	for(int64_t i = 0; i < 16; i++) { \
		_s += sprintf(_s, " %02x", (uint8_t)i); \
	} \
	_s += sprintf(_s, "\n"); \
	for(int64_t i = 0; i < ((len) + 15) / 16; i++) { \
		_s += sprintf(_s, "0x%016" PRIx64 ":", (uint64_t)_ptr); \
		for(int64_t j = 0; j < 16; j++) { \
			_s += sprintf(_s, " %02x", (uint8_t)_ptr[j]); \
		} \
		_s += sprintf(_s, "  "); \
		for(int64_t j = 0; j < 16; j++) { \
			_s += sprintf(_s, "%c", isprint(_ptr[j]) ? _ptr[j] : ' '); \
		} \
		_s += sprintf(_s, "\n"); \
		_ptr += 16; \
	} \
	(char const *)_str; \
})
#endif

#else
#define debug(...) {}
#define dbprintf(fmt, ...) {}
// #define dump(ptr, len) ;
#endif



#define note(...) { \
	note_impl(__VA_ARGS__, ""); \
}
#define note_impl(fmt, ...) { \
	ntprintf("[%s: %s(%d)] " fmt "%s\n", __FILE__, __func__, __LINE__, __VA_ARGS__); \
}
#define ntprintf(fmt, ...) { \
	fprintf(stderr, fmt, __VA_ARGS__); \
}


/* trap */
#undef trap
#ifdef DEBUG_TRAP
#define trap(...) { \
	debug("" __VA_ARGS__); \
	*((volatile uint8_t *)NULL); \
}
#else
#define trap() {}
#endif


/* assertion */
#ifdef DEBUG_ASSERT
#define assert_impl(expr, fmt, ...) { \
	if(!(expr)) { \
		debug("assertion failed (%s) " fmt "%s", #expr, __VA_ARGS__); \
		trap(); \
	} \
}
#define assert(expr, ...) { \
	assert_impl(expr, "" __VA_ARGS__, ""); \
}
#else
#define assert(expr, ...) {}
#endif


/* scan memory for memory sanitizer and valgrind */
#ifdef DEBUG_SCAN
#define _scan_memory(_p, _l) ({ \
	volatile size_t _cnt = 0, _len = (_l); \
	uint8_t const *_ptr = (uint8_t const *)(_p); \
	for(size_t _i = 0; _i < _len; _i++) { if(_ptr[_i]) { _cnt++; } } \
})
#else
#define _scan_memory(_p, _l)	;
#endif


#endif /* #ifndef _DEBUG_H_INCLUDED */
/**
 * end of debug.h
 */
