import json, urllib.request, urllib.parse, sys
qs = {
 'co2-seis-opt': 'CO2 storage injection rate optimization induced seismicity',
 'co2-pressure-mgmt-opt': 'CO2 injection schedule optimization fault reactivation pressure management',
 'convex-seis': 'convex optimization induced seismicity injection',
 'scp-reservoir': 'sequential convex programming reservoir simulation injection optimization',
 'rsf-surrogate': 'rate-and-state surrogate seismicity rate model injection control',
 'cmaes-seis': 'stochastic optimization injection schedule rate-and-state earthquake',
 'dieterich-calib-rsf': 'Dieterich seismicity rate model comparison rate-and-state fault simulation injection',
 'edfm-subcell-leakoff': 'embedded discrete fracture growth leak-off sub-grid storage',
 'gas-frac-growth': 'gas injection hydraulic fracture growth simulation embedded',
 'multiphase-coulomb': 'multiphase CO2 brine Coulomb stress induced seismicity nonlinear',
}
out={}
for tag,s in qs.items():
    url='https://api.openalex.org/works?'+urllib.parse.urlencode({'search':s,'filter':'from_publication_date:2022-01-01','sort':'relevance_score:desc','per_page':10,'select':'id,title,publication_date,doi,cited_by_count','mailto':'research@example.org'})
    try: r=json.load(urllib.request.urlopen(url,timeout=30))['results']
    except Exception as e: print(tag,'ERR',e); continue
    out[tag]=r
    print('##',tag)
    for w in r: print('  ',w['publication_date'],w['cited_by_count'],w['title'][:150],w['doi'] or '')
json.dump(out,open('q47.json','w'))
