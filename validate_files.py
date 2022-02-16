import yaml
import sys
from functools import wraps
from ftplib import FTP
from urllib.parse import urlparse, parse_qsl, urlencode, urlunparse, ParseResult
from urllib.request import Request, urlopen


def retry(number):
    def dec(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for i in range(number):
                try:
                    result = func(*args, **kwargs)
                    return result
                except Exception as e:
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
        return get_dropbox_filesize(url_parsed)

    file = urlopen(Request(url, method='HEAD'))
    return int(file.headers.get('Content-Length'))


def get_dropbox_filesize(url_parsed: ParseResult) -> int:
    # https://stackoverflow.com/a/50067550/4377521
    # otherwise you'll get html page, not file
    url_dict = dict(parse_qsl(url_parsed.query))
    url_dict.update({'dl': 1})
    query = urlencode(url_dict)
    url_parsed = url_parsed._replace(query=query)
    url = urlunparse(url_parsed)
    file = urlopen(Request(url, method='HEAD'))
    return int(file.headers.get('X-Dropbox-Content-Length'))


def get_url(url, access_keys):
    if url.startswith('https://dataset1000genomes.blob.core.windows.net'):
        return url + access_keys['1000g_public_key']
    if url.startswith('https://bioinformatics.file.core.windows.net'):
        return url + access_keys['azure_public_key']
    return url


def get_filesize(data, access_keys):
    if 'expand_rule' in data:
        url = data['url']
        key = data['expand_rule']['key']
        values = data['expand_rule']['values']
        return [get_filesize({'url': url.replace(key, str(val))}, access_keys) for val in values]
    else:
        url = get_url(data['url'], access_keys)

    if urlparse(url).scheme == 'ftp':
        return get_ftp_filesize(url)
    return get_http_filesize(url)


def main(yaml_url):
    with open(yaml_url, 'r') as file:
        content = yaml.safe_load(file)
    access_keys = {k: content.get(k) for k in ('azure_public_key', '1000g_public_key')}
    errors = {}
    for k, v in content['reference'].items():
        if 'url' not in v or 'filesize' not in v:
            continue
        filesize = get_filesize(v, access_keys)
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
