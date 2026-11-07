<%namespace name="common" file="common.cu"/>
${common.defs(preamble, float_type)}
% for name, ev, dim in [("disp_fs", common.disp_fs, 3), ("disp_hs", common.disp_hs, 3), ("strain_fs", common.strain_fs, 6), ("strain_hs", common.strain_hs, 6)]:
static void tde_${name}(const Real* o, const Real* t, const Real* s, Real nu, Real* res) {
    Real3 obs = make3(o[0], o[1], o[2]);
    Real3 tri0 = make3(t[0], t[1], t[2]);
    Real3 tri1 = make3(t[3], t[4], t[5]);
    Real3 tri2 = make3(t[6], t[7], t[8]);
    Real3 slip = make3(s[2], s[0], s[1]);
    ${ev("tri")}
    % for d in range(dim):
    res[${d}] = full_out.${["x","y","z","a","b","c"][d]};
    % endfor
}
% endfor
