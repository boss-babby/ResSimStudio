import json, sys, urllib.request, urllib.parse
def q(search, since='2018-01-01', n=15, field='search'):
    url = 'https://api.openalex.org/works?' + urllib.parse.urlencode({field: search, 'filter': f'from_publication_date:{since}',
          'per_page': n, 'select': 'title,publication_date,doi,cited_by_count', 'mailto': 'research@example.org'})
    d = json.load(urllib.request.urlopen(url, timeout=30))
    return d['meta']['count'], d['results']
if __name__ == '__main__':
    out = {}
    for s in sys.argv[1:]:
        c, r = q(s); out[s] = r
        print(f"## {s}  (count {c})")
        for w in r: print(f"  {w['publication_date']}  [{w['cited_by_count']}] {w['title'][:120]}  {w['doi']}")
    json.dump(out, open('lit/last_query.json', 'w'))
