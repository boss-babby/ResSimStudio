/* C33 C API of the ResSim engine (ressim.dll on Windows, libressim.so on Linux) - used by the Python package `ressim`
 * (ctypes, no compiler needed on the user side) and by any language with a C FFI (MATLAB loadlibrary, Julia ccall,
 * R .C, C#, Excel VBA ...).
 *
 * Conventions: functions returning int give >= 0 on success and -1 on error (rs_last_error() has the text); strings
 * returned are owned by the library and stay valid until the next call on the same session; arrays are copied into
 * caller buffers (call with buf = NULL to get the length). Sessions are independent (one thread per session). */
#ifndef RESSIM_CAPI_H
#define RESSIM_CAPI_H

#if defined(_WIN32)
#  if defined(RESSIM_CAPI_BUILD)
#    define RS_API __declspec(dllexport)
#  else
#    define RS_API __declspec(dllimport)
#  endif
#else
#  define RS_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rs_session rs_session;

RS_API const char* rs_version(void);
RS_API const char* rs_last_error(void);                      /* of the calling thread */
RS_API void rs_set_threads(int n);                            /* OpenMP threads of the engine (0 = all cores) */

/* ---- open / edit */
RS_API rs_session* rs_open(const char* path);                 /* ECLIPSE deck, .rsim project or "example:ID"; NULL on error */
RS_API rs_session* rs_open_text(const char* deck_text, const char* base_dir);
RS_API void rs_close(rs_session* s);
RS_API int rs_set_keyword(rs_session* s, const char* key_value); /* "PERMX=300*200", "RESSIM:FAULTSLIP=..." */
RS_API int rs_num_warnings(rs_session* s);
RS_API const char* rs_warning(rs_session* s, int i);
RS_API const char* rs_info_json(rs_session* s);               /* title, units, engine, dims, wells, report steps, start */
RS_API int rs_list_examples(char* buf, int len);               /* newline-separated ids; returns needed length */

/* ---- run */
RS_API int rs_run(rs_session* s);                             /* whole schedule, any engine */
RS_API int rs_begin(rs_session* s);                           /* black oil: setup + equilibration */
RS_API int rs_step(rs_session* s);                            /* one report step: 1 more, 0 finished, -1 error */
RS_API int rs_apply_schedule(rs_session* s, const char* schedule_text); /* SCHEDULE keywords at the current step */
RS_API double rs_time(rs_session* s);                         /* days */
RS_API int rs_report_step(rs_session* s);
RS_API int rs_num_report_steps(rs_session* s);
RS_API const char* rs_run_json(rs_session* s);                /* ok, error, engine, seconds, stats, kpi, log tail */

/* ---- summary vectors (FOPR, FWCT, FPR, WBHP:P1, ...; "TIME" = days) */
RS_API int rs_summary_count(rs_session* s);
RS_API const char* rs_summary_key(rs_session* s, int i);
RS_API const char* rs_summary_unit(rs_session* s, int i);
RS_API int rs_summary_length(rs_session* s);
RS_API int rs_summary_get(rs_session* s, const char* key, double* buf, int len);

/* ---- per-cell arrays */
RS_API int rs_grid_dims(rs_session* s, int* nx, int* ny, int* nz, int* nactive); /* nactive -1 before begin/run */
RS_API int rs_static_array(rs_session* s, const char* name, double* buf, int len); /* global cells: PORO PERMX ... ACTNUM */
RS_API int rs_active_index(rs_session* s, int* buf, int len);  /* active -> global (0-based) */
RS_API int rs_num_snapshots(rs_session* s);
RS_API double rs_snapshot_time(rs_session* s, int snap);
RS_API const char* rs_cell_array_names(rs_session* s, int snap); /* newline separated */
RS_API int rs_cell_array(rs_session* s, int snap, const char* name, float* buf, int len); /* per active cell */

/* ---- output */
RS_API int rs_write_eclipse(rs_session* s, const char* base);  /* BASE.EGRID/INIT/SMSPEC/UNSMRY/UNRST */
RS_API int rs_export_deck(rs_session* s, const char* path);
RS_API int rs_report(rs_session* s, const char* path, const char* options); /* .docx / .pptx / .html by extension */

/* ---- ECLIPSE / OPM binary summaries and comparison */
RS_API rs_session* rs_open_summary(const char* smspec_path);  /* read-only session holding an external summary */
RS_API const char* rs_compare_json(rs_session* ref, rs_session* cmp, const char* keys_csv); /* per-vector deviations */

/* ---- C43 geomodel scripting (Petrel-parity P10): a geomodel project handle with surfaces, wells, polygons and a
 * seismic volume. Surfaces: z is nx*ny, i fastest, NaN undefined. Indices are 0-based; functions adding an object
 * return its index. JSON strings are owned by the handle (valid until the next call on it). */
typedef struct rs_geo rs_geo;
RS_API rs_geo* rs_geo_new(void);
RS_API void rs_geo_free(rs_geo* g);
RS_API int rs_geo_demo(rs_geo* g, const char* dir, int field_units);          /* writes + loads the synthetic demo project */
/* surfaces */
RS_API int rs_geo_load_surface(rs_geo* g, const char* path);                   /* IRAP classic or ZMAP+ */
RS_API int rs_geo_num_surfaces(rs_geo* g);
RS_API const char* rs_geo_surface_json(rs_geo* g, int i);                      /* name nx ny x0 y0 dx dy rot min max defined */
RS_API int rs_geo_surface_values(rs_geo* g, int i, double* buf, int len);
RS_API int rs_geo_set_surface(rs_geo* g, const char* name, int nx, int ny, double x0, double y0, double dx, double dy,
                              double rot_deg, const double* z);                /* adds, or replaces the same name */
RS_API int rs_geo_surface_sample(rs_geo* g, int i, const double* x, const double* y, int n, double* out);
RS_API int rs_geo_surface_calc(rs_geo* g, const char* expr, int onto, const char* name);   /* names of the surfaces, X, Y, functions */
RS_API int rs_geo_grid_points(rs_geo* g, const double* xyz, int n, int lattice_of, double range, const char* name); /* kriging; lattice_of < 0: own lattice */
RS_API int rs_geo_write_surface(rs_geo* g, int i, const char* path);           /* IRAP classic */
/* wells */
RS_API int rs_geo_load_well_heads(rs_geo* g, const char* path);                /* returns the number of wells added */
RS_API int rs_geo_load_deviation(rs_geo* g, const char* well, const char* path);
RS_API int rs_geo_load_tops(rs_geo* g, const char* path);                      /* returns the number of tops attached */
RS_API int rs_geo_load_las(rs_geo* g, const char* well, const char* path);     /* well NULL/"": name from the LAS header */
RS_API int rs_geo_num_wells(rs_geo* g);
RS_API const char* rs_geo_well_json(rs_geo* g, int i);                         /* name x y kb td, tops [name md x y z], logs [mnemonics] */
RS_API int rs_geo_well_log(rs_geo* g, int i, const char* mnem, double* md, double* v, int len);
RS_API int rs_geo_well_position(rs_geo* g, int i, const double* md, int n, double* xyz);   /* x y tvd per md */
/* polygons */
RS_API int rs_geo_load_polygons(rs_geo* g, const char* path);                  /* returns the number added */
RS_API int rs_geo_add_polygon(rs_geo* g, const char* name, int closed, const double* xy, int n);
RS_API int rs_geo_num_polygons(rs_geo* g);
RS_API const char* rs_geo_polygon_json(rs_geo* g, int i);                      /* name closed area length points */
RS_API int rs_geo_write_polygons(rs_geo* g, const char* path);
RS_API int rs_geo_clip_surface(rs_geo* g, int surface, int polygon, int keep_inside);   /* in place */
/* layer cake: n interfaces (time surfaces by index, well-top names comma separated, per-layer mode 0 mean / 1 map / 2 user
 * and user velocity); depth surfaces are added (named after the tops); returns JSON vmean nwells rms_before rms_after */
RS_API const char* rs_geo_layer_cake(rs_geo* g, const int* time_surfaces, const char* tops_csv, const int* mode, const double* v_user,
                                     int n, double datum, int residual);
/* seismic volume */
RS_API int rs_geo_load_segy(rs_geo* g, const char* path);
RS_API const char* rs_geo_seismic_json(rs_geo* g);                             /* name ni nj ns t0 dt */
RS_API int rs_geo_seismic_trace(rs_geo* g, int i, int j, float* buf, int len);
RS_API int rs_geo_seismic_calc(rs_geo* g, const char* expr, int window);       /* AMP, attributes, T, I, J; replaces the volume */
RS_API int rs_geo_seismic_crop(rs_geo* g, int i0, int i1, int j0, int j1, int s0, int s1, int di, int dj, int ds);
RS_API int rs_geo_write_segy(rs_geo* g, const char* path);

#ifdef __cplusplus
}
#endif
#endif
