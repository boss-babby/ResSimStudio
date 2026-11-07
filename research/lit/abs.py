import json, sys, urllib.request
for doi in sys.argv[1:]:
    try:
        d = json.load(urllib.request.urlopen(f'https://api.openalex.org/works/doi:{doi}?mailto=research@example.org', timeout=30))
        inv = d.get('abstract_inverted_index') or {}
        pos = sorted((p, w) for w, ps in inv.items() for p in ps)
        print(f"### {d['title']} ({d['publication_date']}) {doi}\n{' '.join(w for _, w in pos)[:1600]}\n")
    except Exception as e: print('ERR', doi, e)
