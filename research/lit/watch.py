"""Recent-paper watch: runs saved OpenAlex queries and writes lit/digest_<date>.md (new items only)."""
import json, os, datetime, urllib.request, urllib.parse
Q = {  # topic -> query
 'seismicity-control': 'induced seismicity injection optimization control',
 'rate-state-injection': 'rate-and-state seismicity fluid injection',
 'post-injection': 'post-injection seismicity shut-in strategy',
 'hf-ilsa-edfm': 'hydraulic fracture propagation tip asymptote embedded discrete fracture',
 'phase-field-hf': 'phase-field hydraulic fracture poroelastic',
 'egs-thm': 'enhanced geothermal system thermo-hydro-mechanical simulation',
 'fault-contact-efem': 'embedded finite element fault contact poroelastic reservoir',
 'reservoir-solver': 'nonlinear solver reservoir simulation fully implicit acceleration',
 'adjoint-geomech': 'adjoint optimization geomechanics injection',
 'rsf-schedule': 'rate-and-state fault injection rate schedule event count',
 'hf-leakoff-grid': 'hydraulic fracture growth leak-off grid convergence embedded fracture',
}
here = os.path.dirname(os.path.abspath(__file__)); seenf = os.path.join(here, 'seen.json')
seen = set(json.load(open(seenf))) if os.path.exists(seenf) else set()
since = (datetime.date.today()-datetime.timedelta(days=120)).isoformat(); lines = []
for tag, s in Q.items():
    url = 'https://api.openalex.org/works?'+urllib.parse.urlencode({'search': s, 'filter': f'from_publication_date:{since}',
          'sort': 'relevance_score:desc', 'per_page': 12, 'select': 'id,title,publication_date,doi', 'mailto': 'research@example.org'})
    try: res = json.load(urllib.request.urlopen(url, timeout=30))['results']
    except Exception as e: lines.append(f'## {tag}: ERROR {e}'); continue
    new = [w for w in res if w['id'] not in seen]
    lines.append(f'## {tag} ({len(new)} new)')
    for w in new: lines.append(f"- {w['publication_date']} {w['title']} {w['doi'] or w['id']}"); seen.add(w['id'])
fn = os.path.join(here, f'digest_{datetime.date.today().isoformat()}.md')
open(fn, 'w').write(f'# Literature digest {datetime.date.today()} (since {since})\n\n'+'\n'.join(lines)+'\n')
json.dump(sorted(seen), open(seenf, 'w')); print('wrote', fn)
