#!/usr/bin/env python

import json
from glob import glob

if __name__=='__main__':


    output = open('../../data/nyt/articles.tsv', 'w')
    for filename in glob('../../data/nyt/json/mar*/*.txt'):
        try:
            article = json.load(open(filename, 'r'))
            url = article.keys()[0]
            article = article[url]

            id = article['data']['cms']['article']['id']['value']
            section = article['section']
            title = article['data']['cms']['article']['headline']
            title = title.encode('ascii','replace')

            output.write('%s\t%s\t%s\t%s\n' % (id, section, title, url))

        except ValueError:
            print "error parsing", filename

    output.close()