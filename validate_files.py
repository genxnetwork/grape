import yaml
import sys
from functools import wraps
from ftplib import FTP
from urllib.parse import urlparse, parse_qsl, urlencode, urlunparse
from urllib.request import Request, urlopen


def retry(number):
    def dec(func):
        @wraps
        def wrapper(*args, **kwargs):
            for i in range(number):
                try:
                    result = func(*args, **kwargs)
                    return result
                except Exception:
                    pass
        return wrapper
    return dec


@retry(number=3)
def get_ftp_filesize(url):
    with FTP(urlparse(url).hostname) as ftp:
        ftp.login()
        ftp.sendcmd("TYPE i")  # convert to binary mode
        return ftp.size(urlparse(url).path)


@retry(number=3)
def get_http_filesize(url):
    url_parsed = urlparse(url)
    if 'dropbox' in url_parsed.hostname:
        # https://stackoverflow.com/a/50067550/4377521
        # otherwise you'll get html page, not file
        url_dict = dict(parse_qsl(url_parsed.query))
        url_dict.update({'dl': 1})
        query = urlencode(url_dict)
        url_parsed = url_parsed._replace(query=query)
        url = urlunparse(url_parsed)

    file = urlopen(Request(url, method='HEAD'))
    return file.headers.get('Content-Length')


def get_filesize(data):
    url = data['url']
    if 'expand_rule' in data:
        key = data['expand_rule']['key']
        values = data['expand_rule']['values']
        return [get_filesize({'url': url.replace(key, str(val))}) for val in values]

    if urlparse(url).scheme == 'ftp':
        return get_ftp_filesize(url)
    return get_http_filesize(url)


def main(yaml_url):
    with open(yaml_url, 'r') as file:
        content = yaml.safe_load(file)

    errors = {}
    for k, v in content['reference'].items():
        if 'url' not in v:
            continue
        filesize = get_filesize(v)
        if str(filesize) != str(v['filesize']):
            errors[k] = f'Expected {v["filesize"]}, got {filesize}'

    if errors:
        print(errors)
        sys.exit(2)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Incorrect number of argument, expected 2, got {len(sys.argv)}')
        sys.exit(1)
    main(sys.argv[1])
